import gzip
from concurrent.futures import ThreadPoolExecutor
import logging

logger = logging.getLogger(__name__)


def _open_fastq(file_path: str):
    """Open a FASTQ file, handling gzip compression if needed."""
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')


def read_fastq_in_chunks(file_path: str, chunk_size: int = 1000):
    """
    Generator that yields FASTQ records in chunks.
    Each record consists of 4 lines; this function yields a list of records,
    where each record is a tuple of 4 lines.
    """
    with _open_fastq(file_path) as f:
        chunk = []
        record = []
        for line in f:
            record.append(line.strip())
            if len(record) == 4:
                chunk.append(tuple(record))
                record = []
                if len(chunk) >= chunk_size:
                    logger.debug("Yielding a chunk of %d records from %s", len(chunk), file_path)
                    yield chunk
                    chunk = []
        if chunk:
            logger.debug("Yielding final chunk of %d records from %s", len(chunk), file_path)
            yield chunk


def load_paired_fastq_multithreaded(file_path_r1: str, file_path_r2: str,
                                    num_threads: int = 4, chunk_size: int = 1000):
    """
    Load paired FASTQ files (R1 and R2) using multi-threading.
    Returns a list of paired records where each pair is a tuple:
    (record_R1, record_R2).
    """
    def process_chunk(pair_chunk):
        chunk_r1, chunk_r2 = pair_chunk
        paired = list(zip(chunk_r1, chunk_r2))
        logger.debug("Processed a paired chunk with %d records", len(paired))
        return paired

    paired_reads = []
    # Generators for both FASTQ files.
    r1_chunks = read_fastq_in_chunks(file_path_r1, chunk_size)
    r2_chunks = read_fastq_in_chunks(file_path_r2, chunk_size)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for chunk_r1, chunk_r2 in zip(r1_chunks, r2_chunks):
            futures.append(executor.submit(process_chunk, (chunk_r1, chunk_r2)))
        for future in futures:
            paired_reads.extend(future.result())

    logger.info("Loaded %d paired reads from %s and %s",
                len(paired_reads), file_path_r1, file_path_r2)
    return paired_reads