# csHTO

A bioinformatics tool for demultiplexing single-cell libraries with sample hashtags.

## Features

- **Configuration Validation:** Checks input YAML file integrity (e.g. required sections, existence of R1/R2 FASTQ files).
- **Multi-threaded FASTQ Processing:** Efficiently load large FASTQ files using multi-threading.
- **HTO Demultiplexing & UMI Deduplication:** Extracts HTO sequences, UMIs, and cell barcodes; categorizes reads by sample; deduplicates UMIs.
- **Cell Barcode Extraction & Statistics:** Outputs cell barcode files per sample and detailed statistics.
- **Demultiplexing:** Applies extracted cell barcode filters to additional libraries to produce filtered FASTQ files.
- **Verbose Mode:** Optionally prints detailed step-by-step progress and statistics.

## Installation

1. Clone the repository.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

Run the tool with:

```bash
python -m src.main --input config/sample_config.yaml --output results --threads 4 --verbose
```

## Testing

Run unit tests using:

```python
pytest tests/
```