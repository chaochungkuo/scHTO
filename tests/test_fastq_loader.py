import os
from src.fastq_loader import read_fastq_in_chunks


def test_read_fastq_in_chunks(tmp_path):
    # Create a small FASTQ file.
    fastq_content = "\n".join([
        "@read1",
        "ACGTACGTACGT",
        "+",
        "FFFFFFFFFFFF",
        "@read2",
        "TGCACTGCACTG",
        "+",
        "FFFFFFFFFFFF"
    ])
    file_path = tmp_path / "test.fastq"
    file_path.write_text(fastq_content)
    chunks = list(read_fastq_in_chunks(str(file_path), chunk_size=1))
    print(chunks)
    # Should yield one chunk with 2 records.
    assert len(chunks) == 2
    assert len(chunks[0]) == 1
