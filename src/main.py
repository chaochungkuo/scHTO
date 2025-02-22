import argparse
import os
import logging
from src.config_validator import load_config, ConfigValidationError
from src.fastq_loader import load_paired_fastq_multithreaded
from src.demultiplexer import extract_fields, deduplicate_umis, categorize_reads_by_hto
from src.output_writer import write_cell_barcodes, write_statistics, write_demultiplexed_fastq


def main():
    parser = argparse.ArgumentParser(
        description="Demultiplex single-cell libraries with HTOs."
    )
    parser.add_argument("--input", required=True, help="Path to the configuration YAML file.")
    parser.add_argument("--output", required=True, help="Output directory for results.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    # Set up logging.
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)
    logger.info("Starting cshto with configuration file: %s", args.input)

    # Load and validate configuration.
    try:
        config = load_config(args.input)
        logger.debug("Loaded configuration: %s", config)
    except ConfigValidationError as e:
        logger.error("Configuration Error: %s", e)
        return

    os.makedirs(args.output, exist_ok=True)

    # --- Step 1: Process HTO libraries ---
    # Group FASTQ paths for HTO libraries (expects paired R1 and R2 entries by htolib_name).
    hto_files = {}
    for entry in config["libraries_with_HTOs"]:
        lib = entry["htolib_name"]
        hto_files.setdefault(lib, {})[entry["R12"]] = entry["path"]

    all_extracted_reads = []
    for lib, files in hto_files.items():
        if "R1" not in files or "R2" not in files:
            logger.error("Missing R1 or R2 for library %s", lib)
            continue
        logger.info("Processing HTO library for %s: R1=%s, R2=%s", lib, files["R1"], files["R2"])
        paired_reads = load_paired_fastq_multithreaded(
            files["R1"], files["R2"], num_threads=args.threads
        )
        for pair in paired_reads:
            extracted = extract_fields(pair, config["positions"])
            all_extracted_reads.append(extracted)
        logger.info("Finished processing library %s: %d reads processed", lib, len(paired_reads))

    # --- Step 2: Categorize reads by HTO and perform UMI deduplication ---
    categorized_reads = categorize_reads_by_hto(all_extracted_reads, config["HTO_sequences"])
    sample_cell_barcodes = {}
    dedup_stats = []
    for sample, reads in categorized_reads.items():
        deduped = deduplicate_umis(reads, config["cutoffs"]["min_umi"])
        sample_cell_barcodes[sample] = list(deduped.keys())
        total_umis = sum(len(umis) for umis in deduped.values())
        expected_cells = next(
            (item["estimate_number"] for item in config["expected_cell_number"]
             if item["sample_name"] == sample),
            "NA"
        )
        dedup_stats.append({
            "sample": sample,
            "cell_barcodes": len(deduped),
            "total_unique_umis": total_umis,
            "expected_cells": expected_cells
        })
    logger.info("UMI deduplication complete. Statistics: %s", dedup_stats)

    # --- Step 3: Save cell barcodes and statistics ---
    for sample, cell_barcodes in sample_cell_barcodes.items():
        write_cell_barcodes(args.output, sample, cell_barcodes)
    write_statistics(args.output, dedup_stats)

    # --- Step 4: Demultiplex libraries based on extracted cell barcodes ---
    demux_files = {}
    for entry in config["libraries_to_be_demultiplexed"]:
        lib = entry["htolib_name"]
        demux_files.setdefault(lib, {})[entry["R12"]] = entry["path"]

    for lib, files in demux_files.items():
        if "R1" not in files or "R2" not in files:
            logger.error("Missing R1 or R2 for demultiplexing library %s", lib)
            continue
        logger.info("Processing demultiplexing for library %s: R1=%s, R2=%s", lib, files["R1"], files["R2"])
        paired_reads = load_paired_fastq_multithreaded(
            files["R1"], files["R2"], num_threads=args.threads
        )
        for sample, valid_barcodes in sample_cell_barcodes.items():
            filtered_r1 = []
            filtered_r2 = []
            for pair in paired_reads:
                fields = extract_fields(pair, config["positions"])
                if fields["cell_barcode"] in valid_barcodes:
                    filtered_r1.append(pair[0])
                    filtered_r2.append(pair[1])
            sample_out_dir = os.path.join(args.output, f"demux_{sample}")
            write_demultiplexed_fastq(sample_out_dir, sample, filtered_r1, filtered_r2)
            logger.info("Sample %s: kept %d read pairs after filtering", sample, len(filtered_r1))

    logger.info("Processing complete. Results are saved in: %s", args.output)


if __name__ == "__main__":
    main()
