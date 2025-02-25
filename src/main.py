import argparse
import os
import logging
from src.config_validator import load_config, ConfigValidationError
from src.fastq_loader import load_fastq
from src.demultiplexer import categorize_reads_by_hto, deduplicate_umi, filter_cellbarcodes,split_GEX_fastqs
from collections import OrderedDict


def main():
    parser = argparse.ArgumentParser(
        description="Demultiplex single-cell libraries with HTOs."
    )
    parser.add_argument("--input", required=True, help="Path to the configuration YAML file.")
    parser.add_argument("--output", required=True, help="Output directory for results.")
    parser.add_argument("--threads", type=int, default=6, help="Number of threads to use.")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Chunk size for processing FASTQ files.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    os.environ["MODIN_CPUS"] = str(args.threads)
    # Set up logging.
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt='%Y-%m-%d %H:%M:%S'
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

    for lib, files in hto_files.items():
        samples = [hto["sample_name"] for hto in config["HTO_sequences"]
                   if hto["htolib_name"] == lib]
        statistics = OrderedDict()
        statistics["Library name"] = lib
        if "R1" not in files or "R2" not in files:
            logger.error("Missing R1 or R2 for library %s", lib)
            continue
        statistics["HTO R1 FASTQ"] = files["R1"]
        statistics["HTO R2 FASTQ"] = files["R2"]
        
        # Step 1: Load the libraries with HTOs
        logger.info("Processing HTO library for %s: R1=%s, R2=%s", lib, files["R1"], files["R2"])
        statistics = load_fastq(libname=lib, fastq_R1=files["R1"], fastq_R2=files["R2"],
                   config=config, thread=args.threads,
                   chunk_size=args.chunk_size, output=args.output,
                   statistics=statistics)
        logger.info("Finished processing library %s", lib)
        # Step 2: Categorize read pairs by HTO into samples
        categorized_indices = categorize_reads_by_hto(
            libname=lib, config=config, thread=args.threads,
            chunk_size=args.chunk_size, output=args.output)
        statistics["Valid HTOs"] = sum([len(ind) for ind in categorized_indices.values()])
        for sample in samples:
            statistics[f"{sample} HTOs"] = len(categorized_indices[sample])
        # Step 3: Deduplicate UMIs and categorize cell barcodes
        unique_df = deduplicate_umi(
            categorized_indices=categorized_indices,
            libname=lib,
            thread=args.threads,
            chunk_size=args.chunk_size, 
            output=args.output)
        statistics["Unique barcodes and UMIs"] = unique_df["umi_count"].sum()
        for sample in samples:
            statistics[f"Unique barcodes and UMIs of {sample}"] = unique_df["umi_count"].iloc[unique_df["sample"] == sample].sum()

        # Step 4: Filter cell barcodes from HTO according to the expected cell numbers
        filtered_df = filter_cellbarcodes(unique_df=unique_df,
            config=config, libname=lib,output=args.output)
        statistics["Filtered barcodes"] = filtered_df.shape[0]
        for sample in samples:
            statistics[f"Filtered barcodes of {sample}"] = (filtered_df["sample"] == sample).sum()

        # Step 5: Load GEX library and split the read pairs into samples according to the identified barcodes
        statistics = split_GEX_fastqs(libname=lib, config=config,
                         filtered_df=filtered_df, output=args.output,chunk_size=args.chunk_size,
                         statistics=statistics)
        save_statistics(statistics, output=args.output, libname=lib)
    logger.info("Processing complete. Results are saved in: %s", args.output)


def save_statistics(statistics, output, libname):
    with open(os.path.join(output, f"{libname}_statistics.csv"), "w") as f:
        for k, v in statistics.items():
            print(",".join([k,str(v)]), file=f)
    
if __name__ == "__main__":
    main()
