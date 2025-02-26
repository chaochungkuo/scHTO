import logging
from src.fastq_loader import read_fastq_in_chunks
import numpy as np
import modin.pandas as pd
import os
import concurrent.futures
import gzip
from src.fastq_chunk_processor import process_chunk

logger = logging.getLogger(__name__)


def categorize_reads_by_hto(libname, config, thread,
                            chunk_size, output):
    """
    Categorize reads into samples based on HTO sequences.
    """
    hto_seq = {}
    for hto in config["HTO_sequences"]:
        if hto["htolib_name"] == libname:
            hto_seq[hto["sample_name"]] = hto["hto_sequence"]
    hto_parquet = os.path.join(output, f"{libname}_hto.parquet")
    # Load Parquet file into a DataFrame
    hto_df = pd.read_parquet(hto_parquet)
    logger.info("Load HTO parquet: %s", hto_parquet)
    categorized_indices = {sample_name: [] for sample_name in hto_seq.keys()}

    def process_sample(sample_name, hto_sequence):
        indices = hto_df.index[hto_df['extracted'] == hto_sequence].tolist()
        return sample_name, indices

    with concurrent.futures.ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {executor.submit(process_sample, sample_name, hto_sequence): sample_name for sample_name, hto_sequence in hto_seq.items()}
        for future in concurrent.futures.as_completed(futures):
            sample_name, indices = future.result()
            categorized_indices[sample_name].extend(indices)
            logger.info("Sample %s: %d reads categorized", sample_name, len(indices))
    
    return categorized_indices


def deduplicate_umi(categorized_indices, libname, thread, chunk_size, output):
    bc_parquet = os.path.join(output, f"{libname}_cellbarcodes.parquet")
    umi_parquet = os.path.join(output, f"{libname}_umi.parquet")
    bcumi_df = pd.concat([pd.read_parquet(bc_parquet)['extracted'],
                          pd.read_parquet(umi_parquet)['extracted']], axis=1)
    bcumi_df.columns = ['cell_barcode', 'umi']
    logger.info("Loaded cell barcode and UMI parquet files for %s", libname)
    
    # Add sample column
    sample_arr = np.empty(len(bcumi_df), dtype=object)
    for sample, indices in categorized_indices.items():
        sample_arr[indices] = sample
    bcumi_df['sample'] = sample_arr
    logger.info("Assigned reads to samples for %s", libname)
    # Drop duplicate rows and keep only the unique ones
    # Convert columns to categorical
    bcumi_df['cell_barcode'] = bcumi_df['cell_barcode'].astype('category')
    bcumi_df['umi'] = bcumi_df['umi'].astype('category')
    bcumi_df['sample'] = bcumi_df['sample'].astype('category')
    bcumi_df = bcumi_df.dropna(subset=['sample'])
    unique_df = bcumi_df.drop_duplicates(subset=['cell_barcode', 'umi', 'sample'])
    logger.info("Deduplicated UMIs for %s", libname)
    # Add umi_count for counting the number of unique UMIs per cell barcode per sample
    unique_df = unique_df.groupby(['sample', 'cell_barcode'], observed=True).size().reset_index(name='umi_count')
    logger.info("Calculated the frequency of unique UMIs for %s", libname)
    # Rank rows by umi_count and place the top ones at the top
    unique_df = unique_df.sort_values(by='umi_count', ascending=False)
    logger.info("Ranked cell barcodes by UMI counts for %s", libname)
    
    # Save the result to a file
    unique_df.to_csv(os.path.join(output, f"{libname}_umi_counts.csv"),
                     index=False)
    logger.info("Saved processed %d cell barcodes and UMIs for %s",
                unique_df.shape[0], libname)
    return unique_df

def filter_cellbarcodes(unique_df, config, libname, output):
    cell_numbers = {d["sample_name"]: d["estimate_number"]
                    for d in config["expected_cell_number"]
                    if d["htolib_name"]==libname}
    filtered_df = unique_df.groupby('sample', group_keys=False).apply(
        lambda x: x.nlargest(cell_numbers[x.name], 'umi_count')
    )
    filtered_df.to_csv(os.path.join(output,
        f"{libname}_filtered_cellbarcodes.csv"), index=False)
    logger.info("Saved filtered %d cell barcodes and UMIs for %s",
                filtered_df.shape[0], libname)
    return filtered_df

def split_GEX_fastqs(libname, config, filtered_df, output, chunk_size, statistics):
    # Get dictionary for barcode to sample
    barcode_to_sample = dict(zip(filtered_df['cell_barcode'],
                                 filtered_df['sample']))
    neighbor_dict = build_neighbor_dict(barcode_to_sample)
    # Get cell barcode positions in GEX
    positions = config["positions"]
    cb_start = int(positions['cell_barcode_start']) - 1
    cb_end = int(positions['cell_barcode_end'])
    cb_R12 = positions["cell_barcode_R12"]
    
    samples = [hto["sample_name"] for hto in config["HTO_sequences"]
               if hto["htolib_name"] == libname]
    # Get GEX FASTQ paths
    gex_fastqs = {d["R12"]: d["path"]
                  for d in config["libraries_to_be_demultiplexed"]
                  if d["htolib_name"] == libname}
    
    # Choose the appropriate file open function (gzip or plain text)
    open_func = gzip.open if gex_fastqs["R1"].endswith('.gz') else open

    # Pre-open output files for each sample.
    out_files_R1 = {}
    out_files_R2 = {}
    for sample in samples:
        out_path_R1 = os.path.join(output, f"{libname}_{sample}_R1.fastq.gz")
        out_path_R2 = os.path.join(output, f"{libname}_{sample}_R2.fastq.gz")
        out_files_R1[sample] = gzip.open(out_path_R1, 'wt')
        out_files_R2[sample] = gzip.open(out_path_R2, 'wt')
    # Initialize counters for statistics.
    total_read_pairs = 0
    sample_counts = {sample: 0 for sample in out_files_R1.keys()}
    # Process FASTQ files in lockstep, chunk by chunk.
    with open_func(gex_fastqs["R1"], 'rt') as f1, open_func(gex_fastqs["R2"], 'rt') as f2:
        for chunk1, chunk2 in zip(read_fastq_in_chunks(f1, chunk_size*10),
                                  read_fastq_in_chunks(f2, chunk_size*10)):
            logger.info("Processing a chunk of %d records for %s", len(chunk1), libname)
            local_total, local_counts = process_chunk(chunk1, chunk2,
                                                    barcode_to_sample,
                                                    neighbor_dict,
                                                    cb_start, cb_end, cb_R12,
                                                    out_files_R1, out_files_R2)
            total_read_pairs += local_total
            for sample, count in local_counts.items():
                sample_counts[sample] = sample_counts.get(sample, 0) + count

    # Close all output files.
    for fh in out_files_R1.values():
        fh.close()
    for fh in out_files_R2.values():
        fh.close()
        
    statistics["GEX R1"] = gex_fastqs["R1"]
    statistics["GEX R2"] = gex_fastqs["R2"]
    statistics["GEX total read pairs"] = total_read_pairs
    for sample in set(samples):
        statistics[f"GEX filtered read pairs of {sample}"] = sample_counts[sample]
        
    return statistics

def generate_neighbors(barcode, alphabet="ACGT"):
    neighbors = set()
    barcode = list(barcode)
    for i in range(len(barcode)):
        original = barcode[i]
        for char in alphabet:
            if char != original:
                barcode[i] = char
                neighbors.add("".join(barcode))
        barcode[i] = original  # restore original
    return neighbors

def build_neighbor_dict(barcode_to_sample, alphabet="ACGT"):
    neighbor_dict = {}
    for barcode, sample in barcode_to_sample.items():
        for neighbor in generate_neighbors(barcode, alphabet):
            # Only add if this neighbor isn't already in the exact dictionary,
            # or mark as ambiguous if needed.
            if neighbor not in barcode_to_sample:
                if neighbor in neighbor_dict and neighbor_dict[neighbor] != sample:
                    neighbor_dict[neighbor] = None  # mark ambiguous
                else:
                    neighbor_dict[neighbor] = sample
    return neighbor_dict