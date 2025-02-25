import logging
import os
import modin.pandas as pd
from fastqio import FASTQReader


logger = logging.getLogger(__name__)


def load_fastq(libname, fastq_R1, fastq_R2, config, thread, chunk_size, output, statistics):
    """
    Load R1 and R2 FASTQ files and extract cell barcodes, UMI and HTOs accordingly.
    """
    positions = config["positions"]
    # Cell barcodes
    cb_start = int(positions['cell_barcode_start']) - 1
    cb_end = int(positions['cell_barcode_end'])
    prefix = os.path.join(output, f"{libname}_cellbarcodes")
    if positions["cell_barcode_R12"] == "R1":
        fastq = FASTQReader(fastq_R1, thread=thread, chunk_size=chunk_size)
    else:
        fastq = FASTQReader(fastq_R2, thread=thread, chunk_size=chunk_size)
    fastq.extract(start=cb_start, end=cb_end,
                  save_parquet=True, parquet_prefix=prefix)
    # UMI
    umi_start = int(positions['umi_start']) - 1
    umi_end = int(positions['umi_end'])
    prefix = os.path.join(output, f"{libname}_umi")
    if positions["umi_R12"] == "R1":
        fastq = FASTQReader(fastq_R1, thread=thread, chunk_size=chunk_size)
    else:
        fastq = FASTQReader(fastq_R2, thread=thread, chunk_size=chunk_size)
    fastq.extract(start=umi_start, end=umi_end,
                  save_parquet=True, parquet_prefix=prefix)
    # HTO
    hto_start = int(positions['hto_start']) - 1
    hto_end = int(positions['hto_end'])
    prefix = os.path.join(output, f"{libname}_hto")
    if positions["hto_R12"] == "R1":
        fastq = FASTQReader(fastq_R1, thread=thread, chunk_size=chunk_size)
    else:
        fastq = FASTQReader(fastq_R2, thread=thread, chunk_size=chunk_size)
    fastq.extract(start=hto_start, end=hto_end,
                  save_parquet=True, parquet_prefix=prefix)
    statistics["Total read pair"] = fastq.count_reads()
    
    return statistics


def extract_information(pair, positions, hto_sequences):
    """
    Extracts and encodes hashtag, cell barcode, and UMI from a read pair.
    """
    r1, r2 = pair
    cb_start = int(positions['cell_barcode_start']) - 1
    cb_end = int(positions['cell_barcode_end'])
    if positions["cell_barcode_R12"] == "R1":
        cell_barcode = r1[cb_start:cb_end]
    else:
        cell_barcode = r2[cb_start:cb_end]

    umi_start = int(positions['umi_start']) - 1
    umi_end = int(positions['umi_end'])
    if positions["umi_R12"] == "R1":
        umi = r1[umi_start:umi_end]
    else:
        umi = r2[umi_start:umi_end]

    hashtag_start = int(positions['hto_start']) - 1
    hashtag_end = int(positions['hto_end'])
    if positions["hto_R12"] == "R1":
        hashtag = r1[hashtag_start:hashtag_end]
    else:
        hashtag = r2[hashtag_start:hashtag_end]

    return {
        'cell_barcode': encode_seq(cell_barcode),
        'umi': encode_seq(umi),
        'hashtag': encode_seq(hashtag)
    }

def save_to_parquet(data, output_dir, chunk_index):
    """
    Save extracted information to a Parquet file.
    """
    df = pd.DataFrame(data)
    output_path = os.path.join(output_dir, f"extracted_reads_chunk_{chunk_index}.parquet")
    df.to_parquet(output_path)
    logger.info("Saved chunk %d to %s", chunk_index, output_path)


def read_fastq_in_chunks(file_handle, chunk_size=10000):
    """
    Generator that yields a list of FASTQ records.
    Each record is a tuple: (header, sequence, plus, quality).
    """
    records = []
    while True:
        # Read one record (4 lines)
        record = [file_handle.readline().rstrip() for _ in range(4)]
        if not record[0]:
            # End of file; yield any remaining records.
            if records:
                yield records
            break
        records.append(tuple(record))
        if len(records) >= chunk_size:
            yield records
            records = []