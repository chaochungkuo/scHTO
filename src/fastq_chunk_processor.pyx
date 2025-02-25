# fastq_chunk_processor.pyx
# cython: boundscheck=False, wraparound=False, nonecheck=False, language_level=3

def process_chunk(list chunk1, list chunk2,
                  dict barcode_to_sample,
                  dict neighbor_dict,
                  int cb_start, int cb_end, str cb_R12,
                  dict out_files_R1, dict out_files_R2):
    """
    Process a chunk of FASTQ records, writing paired records to sample files,
    and accumulate counts for statistics.

    Parameters:
      chunk1: List of tuples (header, seq, plus, qual) from R1.
      chunk2: List of tuples (header, seq, plus, qual) from R2.
      barcode_to_sample: Dictionary mapping cell barcode to sample.
      cb_start: Start index for barcode extraction (0-indexed).
      cb_end: End index for barcode extraction.
      cb_R12: R1 or R2 where barcode locates.
      out_files_R1: Dictionary mapping sample to output file object for R1.
      out_files_R2: Dictionary mapping sample to output file object for R2.
    
    Returns:
      A tuple (local_total, chunk_sample_counts) where:
        local_total is the number of paired records processed.
        chunk_sample_counts is a dict mapping sample to the number of records written.
    """
    cdef int i, n = len(chunk1)
    cdef tuple rec1, rec2
    cdef str header1, seq1, plus1, qual1
    cdef str header2, seq2, plus2, qual2
    cdef str cell_barcode, sample
    cdef int local_total = 0

    # Initialize a Python dictionary for sample counts
    chunk_sample_counts = {}
    
    for i in range(n):
        local_total += 1
        rec1 = chunk1[i]
        rec2 = chunk2[i]
        # Unpack R1
        header1 = rec1[0]
        seq1    = rec1[1]
        plus1   = rec1[2]
        qual1   = rec1[3]
        # Unpack R2
        header2 = rec2[0]
        seq2    = rec2[1]
        plus2   = rec2[2]
        qual2   = rec2[3]

        # Extract the cell barcode from the R1 sequence
        if cb_R12 == "R1":
            cell_barcode = seq1[cb_start:cb_end]
        else:
            cell_barcode = seq2[cb_start:cb_end]
        sample = get_sample_for_barcode_fast(cell_barcode,
                                             barcode_to_sample,
                                             neighbor_dict)
        if sample is not None:
            # Update the count for this sample
            if sample in chunk_sample_counts:
                chunk_sample_counts[sample] = chunk_sample_counts[sample] + 1
            else:
                chunk_sample_counts[sample] = 1
            # Write the paired records to the respective sample files.
            out_files_R1[sample].write(header1 + "\n" + seq1 + "\n" + plus1 + "\n" + qual1 + "\n")
            out_files_R2[sample].write(header2 + "\n" + seq2 + "\n" + plus2 + "\n" + qual2 + "\n")

    return local_total, chunk_sample_counts

def get_sample_for_barcode_fast(str cell_barcode, dict barcode_to_sample, dict neighbor_dict):
    """
    Return the sample corresponding to cell_barcode using an exact match first,
    then falling back to the neighbor dictionary.
    
    Parameters:
      cell_barcode: The barcode string.
      barcode_to_sample: Dictionary mapping barcodes to samples (exact matches).
      neighbor_dict: Dictionary mapping barcode neighbors (within Hamming distance 1)
                     to samples.
    
    Returns:
      The sample if found, otherwise None.
    """
    cdef object sample = barcode_to_sample.get(cell_barcode)
    if sample is not None:
        return sample
    sample = neighbor_dict.get(cell_barcode)
    if sample is None:
        return None
    return sample