# bcParse.py

**bcParse.py** is a Python script for parsing short-read sequencing data (FASTQ files) and associated run information (XLSX), then outputting formatted Excel workbooks and long-format CSV files. It is designed for flexible, stock-specific workflows in sequencing data analysis, and supports both command-line and graphical user interface (GUI) modes.

## Key Features

- Multiple stock settings for sequencing data (see `config.py` for stock definitions)
- Dual interface: Command-line and Tkinter-based GUI for file selection
- Excel output: Sheets arranged by sequencing index groups as defined in run info
- Dual-index detection: Supports dual-indexed sequencing (enable with `--dualindex`)
- Quality filtering: Set thresholds for mean Phred quality and mask low-quality bases
- Customizable output: Legacy format and matrix filtering options available

## Requirements

- Python 3.x
- Local modules:
  - `config.py`
  - `utils/io.py`
  - `utils/countdata_parser.py`
  - `utils/fastq_parser.py`
  - `utils/gui.py`
- Third-party packages:
  - `pandas`
  - `tqdm`
  - `numpy`
  - `levenshtein`
- Pre-configured `ref` directory containing necessary reference data.

## Installation

**Make a local copy of the base directory `bcParse_v2.53`.**

Ensure all required third-party packages are installed.  
You can do this by running:

```python
python /path/to/bcParse_v2.53/utils/third_party_packages.py
```

Make sure .py permissions are set.

```bash
find /path/to/bcParse_v2.53/ -type f -name "*.py" -exec chmod u+x {} +
```

## Usage

### Command-Line Mode

```python
python /path/to/bcParse_v2.53/bcParse.py --sample_path path/to/your.fastq --runinfo_path path/to/runinfo.xlsx --out_path path/to/output --stock 239M --dualindex
```

### GUI Mode

We recommend this over CLI. 

```python
python /path/to/bcParse_v2.53/bcParse.py --gui
```

*This opens a file dialog for interactive selection of input/output files and parameters.*

## Arguments

| Argument           | Type    | Description                                                                                 |
|--------------------|---------|---------------------------------------------------------------------------------------------|
| --sample_path      | str     | Path to input FASTQ file                                                                    |
| --runinfo_path     | str     | Path to run info file (CSV/XLSX)                                                            |
| --out_path         | str     | Output directory                                                                            |
| --stock            | str     | Stock-specific run mode (see `config.py` for options)                                       |
| --dualindex        | flag    | Enable dual-index detection                                                                 |
| --mean_qual        | int     | Minimum allowed mean Phred quality per read/sequence (default: 30)                          |
| --mismatches       | int     | Number of allowed mismatches in sequence extraction (default: 1)                            |
| --mask             | flag    | Mask low-quality bases                                                                      |
| --mask_qual        | int     | Quality threshold for masking bases (default: 20)                                           |
| --legacy_format    | flag    | Output in legacy format                                                                     |
| --filt_mat_ac      | flag    | Filter AC in sample group matrix                                                            |
| --gui              | flag    | Launch GUI for file selection                                                               |

## Example

```python
python bcParse.py --sample_path data/sample.fastq --runinfo_path data/runinfo.xlsx --out_path results/ --stock 239M --mean_qual 35 --mismatches 2 --dualindex
```


## Workflow Overview

1. **Initialize arguments and settings**  
   - Parse CLI or GUI inputs.
   - Load stock-specific settings from `config.py`.

2. **Stream and parse FASTQ data**  
   - Count barcode reads per index.
   - Apply quality and mismatch filters.

3. **Output results**  
   - Write formatted Excel workbook and CSV.
   - Arrange Excel sheets by index groups from run info.

## Notes

- Required arguments: `--sample_path`, `--runinfo_path`, `--out_path`, `--stock`
- For new stock configurations, edit `config.py`
- For troubleshooting or custom workflows, review utility modules in `utils/`
- To ensure all dependencies are installed, run `python third_party_packages.py` before first use

## Author & Version

- Author: Charles Goodman
- Version: 2.53 (2025/06/04)

Keele Lab, ACVP, FNLCR
