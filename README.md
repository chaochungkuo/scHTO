# csHTO

A bioinformatics tool for demultiplexing single-cell libraries with sample hashtags.

## Main steps:

1. Load the read pairs of HTO libraries.
2. Extract HTO sequences and assign the read pairs into the defined samples. In this step, cell barcodes and UMIs are also assigned to the samples together with HTO sequences.
3. Deduplicate UMIs per cell barcode per sample.
4. Calculate UMI counts per cell barcodes per sample and rank cell barcodes according to UMI counts.
5. Filter the top cell barcodes by UMI count with the expected cell number per sample.
6. Load the GEX libraries and identify the categorized cell barcodes and export separate FASTQ files for each sample.

## Features

- **Configuration Validation:** Checks input YAML file integrity (e.g. required sections, existence of R1/R2 FASTQ files).
- **Multi-threaded FASTQ Processing:** Efficiently load large FASTQ files using multi-threading.
- **HTO Demultiplexing & UMI Deduplication:** Extracts HTO sequences, UMIs, and cell barcodes; categorizes reads by sample; deduplicates UMIs.
- **Cell Barcode Extraction & Statistics:** Outputs cell barcode files per sample and detailed statistics.
- **Demultiplexing:** Applies extracted cell barcode filters to additional libraries to produce filtered FASTQ files.
- **Hamming distance:** Only cell barcode is allowed for hamming distance 1, because cell barcodes are longer and designed to be distinct. Hamming distance is now applied in matching HTO sequences and UMIs.
- **Verbose Mode:** Optionally prints detailed step-by-step progress and statistics.

## Installation

1. Clone the repository.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Please also install [fastqio](https://github.com/chaochungkuo/fastqio) from the source.
4. Install scHTO in the scHTO folder:
   ```bash
   pip install .
   ```

## Usage

Run the tool with:

```bash
schto --input config/sample_config.yaml --output results --threads 4 --verbose
```

Here is the `sample_config.yaml` which allows you to define all the parameters.

```yaml
libraries_with_HTOs:
  - htolib_name: pool1
    R12: R1
    path: Pool1_cell_surface_protein_sample_R1.fastq.gz
  - htolib_name: pool1
    R12: R2
    path: Pool1_cell_surface_protein_sample_R2.fastq.gz

positions:
  cell_barcode_R12: R1
  cell_barcode_start: 1
  cell_barcode_end: 16
  umi_R12: R1
  umi_start: 17
  umi_end: 28
  hto_R12: R2
  hto_start: 1
  hto_end: 15

libraries_to_be_demultiplexed:
  - htolib_name: pool1
    R12: R1
    path: Pool1_GEX_sample_R1.fastq.gz
  - htolib_name: pool1
    R12: R2
    path: Pool1_GEX_sample_R2.fastq.gz

HTO_sequences:
  - htolib_name: pool1
    sample_name: sample1
    hto_sequence: TTGGCCTTTGTATCG
  - htolib_name: pool1
    sample_name: sample2
    hto_sequence: AACGCCAGTATGAAC
  - htolib_name: pool1
    sample_name: sample3
    hto_sequence: GCTTCCGTATATCTG

expected_cell_number:
  - htolib_name: pool1
    sample_name: sample1
    estimate_number: 30000
  - htolib_name: pool1
    sample_name: sample2
    estimate_number: 30000
  - htolib_name: pool1
    sample_name: sample3
    estimate_number: 30000
```

## Statistics

For transparency, scHTO exports the following files:

- `output/pool1_umi_counts.csv` for UMI counts per cell barcodes per sample where you extracted from the library with HTOs.
  ```csv
  sample,cell_barcode,umi_count
  sample1,CAATCCTGTCTGTGTG,198
  sample3,CTGGAACGTCAGTCCA,159
  sample1,GTTCATCGTCAGCGGA,46
  sample1,CTCGAATTCCGGGTTG,45
  sample3,AAGAGGTGTCATGGCC,42
  sample3,GCGTAGTGTGCAGGTG,33
  sample1,CCCAGCTCATCACATC,33
  sample1,TGGATGTAGCCGAGTT,30
  sample1,GGAAATGTCCATCGTC,30
  sample1,GTGTAGCGTCTGTGCA,28
  sample3,AACGCATAGCATCCTT,27
  sample3,GGGCCATAGCCTGGTG,27
  ...
  ```
- `output/pool1_statistics.csv` for all statistics.
  ```csv
  Library name,pool1
  HTO R1 FASTQ,Pool1_cell_surface_protein_sample_R1.fastq.gz
  HTO R2 FASTQ,Pool1_cell_surface_protein_sample_R2.fastq.gz
  Total read pair,25000
  Valid HTOs,19330
  sample1 HTOs,6263
  sample2 HTOs,2271
  sample3 HTOs,10796
  Unique barcodes and UMIs,17016
  Unique barcodes and UMIs of sample1,5492
  Unique barcodes and UMIs of sample2,2012
  Unique barcodes and UMIs of sample3,9512
  Filtered barcodes,14418
  Filtered barcodes of sample1,4333
  Filtered barcodes of sample2,1883
  Filtered barcodes of sample3,8202
  GEX R1,Pool1_GEX_sample_R1.fastq.gz
  GEX R2,Pool1_GEX_sample_R2.fastq.gz
  GEX total read pairs,250000
  GEX filtered read pairs of sample2,56
  GEX filtered read pairs of sample3,467
  GEX filtered read pairs of sample1,112
  ```

## Testing

Run unit tests using:

```python
pytest tests/
```