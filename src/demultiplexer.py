import logging

logger = logging.getLogger(__name__)


def _build_positions_dict(positions_list):
    """
    Convert positions list from config into a dictionary keyed by (R12, name).
    Example: {("R1", "cell_barcode_start"): 1, ...}
    """
    pos_dict = {}
    for pos in positions_list:
        key = (pos["R12"], pos["name"])
        pos_dict[key] = pos["position"]
    return pos_dict


def extract_fields(read_pair, positions_list):
    """
    Extract cell barcode, UMI, and HTO sequence from a paired read.
    
    Args:
        read_pair: Tuple (record_R1, record_R2), where each record is a tuple of 4 lines.
        positions_list: List of position entries from config.
    
    Returns:
        A dictionary with keys: 'cell_barcode', 'UMI', and 'HTO'.
    """
    pos_dict = _build_positions_dict(positions_list)
    r1 = read_pair[0][1]  # Sequence from R1.
    r2 = read_pair[1][1]  # Sequence from R2.

    try:
        cell_barcode = r1[pos_dict[("R1", "cell_barcode_start")] - 1:
                          pos_dict[("R1", "cell_barcode_end")]]
        umi = r1[pos_dict[("R1", "umi_start")] - 1:
                 pos_dict[("R1", "umi_end")]]
        hto = r2[pos_dict[("R2", "hto_start")] - 1:
                 pos_dict[("R2", "hto_end")]]
    except KeyError as e:
        logger.error("Missing position info: %s", e)
        raise

    logger.debug("Extracted cell_barcode: %s, UMI: %s, HTO: %s", cell_barcode, umi, hto)
    return {"cell_barcode": cell_barcode, "UMI": umi, "HTO": hto}


def deduplicate_umis(reads, min_umi: int):
    """
    Deduplicate UMIs for each cell barcode.
    
    Args:
        reads: List of dictionaries, each with keys 'cell_barcode', 'UMI', and 'HTO'.
        min_umi: Minimum number of unique UMIs required to keep a cell barcode.
    
    Returns:
        A dictionary mapping cell barcode to a set of unique UMIs that meet the minimum count.
    """
    deduped = {}
    for read in reads:
        cb = read["cell_barcode"]
        umi = read["UMI"]
        deduped.setdefault(cb, set()).add(umi)
    # Filter cell barcodes that do not meet the minimum UMI count.
    filtered = {cb: umis for cb, umis in deduped.items() if len(umis) >= min_umi}
    logger.info("Deduplicated UMIs: kept %d cell barcodes after filtering with min_umi=%d",
                len(filtered), min_umi)
    return filtered


def categorize_reads_by_hto(reads, hto_sequences):
    """
    Categorize reads into samples based on HTO sequences.
    
    Args:
        reads: List of read dictionaries (with key 'HTO').
        hto_sequences: List of dicts from config with keys 'sample_name' and 'hto_sequence'.
    
    Returns:
        A dictionary mapping sample_name to a list of read dictionaries.
    """
    categorized = {ht["sample_name"]: [] for ht in hto_sequences}
    for read in reads:
        assigned = False
        for ht in hto_sequences:
            if read["HTO"] == ht["hto_sequence"]:
                categorized[ht["sample_name"]].append(read)
                assigned = True
                break
        if not assigned:
            logger.debug("Read with HTO %s not assigned to any sample", read["HTO"])
    for sample, sample_reads in categorized.items():
        logger.info("Sample %s: %d reads categorized", sample, len(sample_reads))
    return categorized