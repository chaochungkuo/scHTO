import os
import csv
import logging

logger = logging.getLogger(__name__)


def write_cell_barcodes(output_dir, sample, cell_barcodes):
    """
    Write a list of cell barcodes to a file for the given sample.
    """
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, f"{sample}_cell_barcodes.txt")
    with open(file_path, 'w') as f:
        for barcode in cell_barcodes:
            f.write(f"{barcode}\n")
    logger.info("Wrote %d cell barcodes for sample %s to %s",
                len(cell_barcodes), sample, file_path)


def write_statistics(output_dir, statistics):
    """
    Write statistics to a CSV file.
    
    Args:
        statistics: List of dictionaries where keys are statistic names.
    """
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, "statistics.csv")
    if not statistics:
        logger.warning("No statistics to write.")
        return
    keys = statistics[0].keys()
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        writer.writerows(statistics)
    logger.info("Wrote statistics to %s", file_path)


def write_demultiplexed_fastq(output_dir, sample, reads_r1, reads_r2):
    """
    Write demultiplexed FASTQ files for a sample.
    
    Args:
        reads_r1: List of FASTQ records (each a tuple of 4 lines) for R1.
        reads_r2: List of FASTQ records for R2.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_r1 = os.path.join(output_dir, f"{sample}_R1.fastq")
    out_r2 = os.path.join(output_dir, f"{sample}_R2.fastq")
    
    with open(out_r1, 'w') as f1, open(out_r2, 'w') as f2:
        for record in reads_r1:
            f1.write("\n".join(record) + "\n")
        for record in reads_r2:
            f2.write("\n".join(record) + "\n")
    logger.info("Wrote demultiplexed FASTQ for sample %s", sample)