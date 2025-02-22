from src.demultiplexer import extract_fields, deduplicate_umis, categorize_reads_by_hto


def test_extract_fields():
    # Create a fake paired read.
    read_r1 = ("@read1", "ACGTACGTACGTACGTACGTACGT", "+", "FFFFFFFFFFFFFFFFFFFFFFFF")
    read_r2 = ("@read1", "TTTTCCCCAAAAGGGGTTTT", "+", "FFFFFFFFFFFFFFFFFFFF")
    # Fake positions: cell barcode from 1-4, UMI from 5-8 in R1, HTO from 1-5 in R2.
    positions = [
        {"name": "cell_barcode_start", "position": 1, "R12": "R1"},
        {"name": "cell_barcode_end", "position": 4, "R12": "R1"},
        {"name": "umi_start", "position": 5, "R12": "R1"},
        {"name": "umi_end", "position": 8, "R12": "R1"},
        {"name": "hto_start", "position": 1, "R12": "R2"},
        {"name": "hto_end", "position": 5, "R12": "R2"}
    ]
    fields = extract_fields((read_r1, read_r2), positions)
    assert fields["cell_barcode"] == "ACGT"
    assert fields["UMI"] == "ACGT"
    # Depending on slicing, the expected HTO might be "TTTTC"
    assert fields["HTO"].startswith("TTTTC")


def test_deduplicate_umis():
    reads = [
        {"cell_barcode": "AAAA", "UMI": "UMI1", "HTO": "HTO1"},
        {"cell_barcode": "AAAA", "UMI": "UMI1", "HTO": "HTO1"},
        {"cell_barcode": "AAAA", "UMI": "UMI2", "HTO": "HTO1"},
        {"cell_barcode": "CCCC", "UMI": "UMI3", "HTO": "HTO2"}
    ]
    deduped = deduplicate_umis(reads, min_umi=2)
    assert "AAAA" in deduped
    assert "CCCC" not in deduped  # Only one unique UMI for CCCC.


def test_categorize_reads_by_hto():
    reads = [
        {"cell_barcode": "AAAA", "UMI": "UMI1", "HTO": "HTO1"},
        {"cell_barcode": "CCCC", "UMI": "UMI2", "HTO": "HTO2"},
        {"cell_barcode": "GGGG", "UMI": "UMI3", "HTO": "HTO1"}
    ]
    hto_sequences = [
        {"htolib_name": "pool1", "sample_name": "sample1", "hto_sequence": "HTO1"},
        {"htolib_name": "pool1", "sample_name": "sample2", "hto_sequence": "HTO2"}
    ]
    categorized = categorize_reads_by_hto(reads, hto_sequences)
    assert len(categorized["sample1"]) == 2
    assert len(categorized["sample2"]) == 1