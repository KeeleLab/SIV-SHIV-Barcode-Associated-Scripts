#! /usr/bin/env python3
# bcParse.py

# %% Imports

# Standard
import argparse
from argparse import RawTextHelpFormatter

# Local
from config import stocks, get_settings
from utils.io import (
    read_csv_to_dict,
    format_refdata,
    process_fastq_stream,
    write_output,
)
from utils.countdata_parser import parse_countdata


# %% CLI args

ver = "2.53 - 2025/06/04"

description = f"""
Name:       bcParse.py
Author:     CAG
Version:    {ver}

Given a short-read fastq file and associated runinfo csv, parses data into a formatted excel workbook
and long-format csv file. 

Key Features:
- Supports multiple stock settings for sequencing data.
- Supports both command-line and graphical user interface (GUI) modes.
- Arranges excel sheets by sequencing index groups defined in runinfo.
- V2 Now includes dual-index detection for M/M2, when --dualindex is specified.

Usage:
Run the script in CLI mode by default or use the --gui flag to launch the Tkinter GUI.

Example CLI command:
python bcParse_v1.0.py --sample_path path/to/your.fastq --runinfo_path path/to/runinfo.xlsx --out_path path/to/wherever --stock 239M --dualindex

Example GUI command:
python bcParse_v1.0.py --gui # will open a dialogue window for moused selection
"""

# Create the parser
parser = argparse.ArgumentParser(
    description=description, formatter_class=RawTextHelpFormatter
)

# File args
file_args = [
    ("--sample_path", "path/to/sample.fq"),
    ("--runinfo_path", "path/to/runinfo.xlsx"),
    ("--out_path", "path/to/output_location"),
]

for arg, help_text in file_args:
    parser.add_argument(arg, type=str, help=help_text)

# Stock arg
parser.add_argument(
    "--stock", choices=stocks, help="Declare stock-specific run mode"
    )

# Set mean q33 threshold - applied to read & extracted seqs
parser.add_argument(
    "--mean_qual", type=int, default=30, help="Minimum allowed mean phred per read/seq extract"
    )

# Set number of allowed mismatches in seq extraction
parser.add_argument(
    "--mismatches", type=int, default=1, help="Number of allowed target mismatches"
    )

# Mask low quality bases
parser.add_argument(
    "--mask", action="store_true", help="Mask low quality bases"
    )

# Set mask q33 threshold - applied when --mask = T
parser.add_argument(
    "--mask_qual", type=int, default=20, help="Mask low quality bases (default: 20)"
    )

# Legacy format
parser.add_argument(
    "--legacy_format", action="store_true", help="Return output in legacy format"
    )

# Filter AC in samp_group matrix
parser.add_argument(
    "--filt_mat_ac", action="store_true", help="Mask low quality bases"
    )

# GUI arg
parser.add_argument(
    "--gui", action="store_true", help="Use GUI for file selection"
    )

# Parse arguments
args = parser.parse_args()

# %% Main process function


def main():
    """
    1. Initialize arguments and run settings. See config.py, utils/io.py, utils/gui.py

    - Required args are sample_path, runinfo_path, out_path, and stock.
    - See config.py for stock options, e.g. "239M"
    - Remaining args have defaults set - see above. Modify if you know what you're doing! 
    """
    args = parser.parse_args()

    # If running with --gui...
    if args.gui:
        from utils.gui import ArgSelectionGUI

        # set variables returned from tkinter object,
        (
            sample_path,
            runinfo_path,
            out_path,
            stock,
            mean_qual,
            mismatches,
            mask,
            mask_qual,
            filt_mat_ac,
            legacy_format,

        ) = ArgSelectionGUI(stocks).run()

    # otherwise use CLI args.
    else:
        sample_path = args.sample_path
        runinfo_path = args.runinfo_path
        out_path = args.out_path
        stock = args.stock
        mean_qual = args.mean_qual
        mismatches = args.mismatches
        mask = args.mask
        mask_qual = args.mask_qual
        filt_mat_ac = args.filt_mat_ac
        legacy_format = args.legacy_format

    if not all([sample_path, runinfo_path, out_path]):
        print(
            "Error: missing args. Call with -h to see full args, or with --gui for graphic interface."
        )
        return

    # Show off...
    print("""
                                                                         
    ,---.                        |         ,---.                         
    |---.,---.,---.,---.,---.,---|,---.    |---',---.,---.,---.,---.,---.
    |   |,---||    |    |   ||   ||---'    |    ,---||    `---.|---'|    
    `---'`---^`    `---'`---'`---'`---'    `    `---^`    `---'`---'`    
    Keele lab, ACVP, FNLCR                                               
    
    """)

    # Get stock-defined run settings - add new stocks and settings in config.py!
    settings = get_settings(stock)
    
    # Add analysis settings from args:
    settings.update({
        'mismatches': mismatches, 
        'mean_qual': mean_qual,
        'mask' : mask, 
        'mask_qual': mask_qual
    })
    # Write to terminal
    print("\n".join(f"{key}: {value}" for key, value in settings.items()) + "\n")

    """
    2. Stream fastq to count barcodes reads per index
    
    - Return {idx:{bc:count}}
    - See utils/io.py and utils/fastq_parser.py
    - mask returns N for base q<mask_qual
    - reads drop when mean phred<mean_qual
    - mismatches = # allowed errors
    """
    res_pfq = process_fastq_stream(
        sample_path=sample_path, 
        settings=settings, 
    )

    print("\n")

    """
    3. Main analysis: parse {idx_seq:{bc_seq:count}} to multiattr obj res_pcd
    
    - Incl. final processed 'df' & class methods for output 
    - See utils/countdata_parser.py 
    """
    # Format reference data for parse_countdata()
    runinfo, p5_refdict, bc_refdict = format_refdata(
        runinfo_path, settings["primer_path"], settings["barcode_path"]
    )

    # adjust for dual index runs:
    if settings["dualindex"]:
        # Initialize p7_refdict for sec.3
        p7_refdict = read_csv_to_dict(settings["p7_path"])
        # Define full index name for sec.3
        runinfo["full_idx"] = runinfo["Barcodes"] + "_" + runinfo["(F Barcode)"]
    else:
        # Default to p5-only settings
        p7_refdict = None
        runinfo["full_idx"] = runinfo["Barcodes"]

    # Run and hold instance as res_pcd:
    res_pcd = parse_countdata(
        countdict=res_pfq,
        runinfo=runinfo,
        bc_refdict=bc_refdict,
        p5_refdict=p5_refdict,
        p7_refdict=p7_refdict,
    )

    print("\n")

    """ 
    4. Main output - see io.py & utils/countdata_parser.py. 
    """
    # Collect remaining args in dict to pass to output
    args_dict = {
        "sample_path": sample_path,
        "runinfo_path": runinfo_path,
        "out_path": out_path,
        "stock": stock,
        "filt_mat_ac": filt_mat_ac,
        "legacy_format": legacy_format,
        "software version": ver,
    }

    write_output(
        runinfo_df=runinfo,
        settings_dict=settings,
        args_dict=args_dict,
        res_pcd=res_pcd,
        out_path=out_path,
        legacy_format=legacy_format,
        filt_mat_ac=filt_mat_ac,
    )

    # In case of -i run, to examine:
    return res_pcd

if __name__ == "__main__":
    result = main()

# %%
""" 
Version history

v2.20:
1. parse_fastq_rec returned as instance instead of dict repr
2. p5/p7_seq and p5/p7_qual are held as independent attrs during parse_fastq_rec init, joined (or not) in process_fastq_stream()
3. bases < q20 are optionally substituted with N, resulting in idx and bc_seqs showing N, which can be leveraged for finer-grain analyses
4. --mask, --filt_mat_ac, --legacy_format added as optional boolean params, default False
5. GUI modified to indicate new params

v2.30
1. removed --dualindex as argument and defined as per-stock setting in config.py (incl updates to gui and args)
2. added 239M/M2_dualindex as stock settings to accomodate 1.
3. BIG: for 239M/M2_dualindex, I extract VPX to VPR as barcode, whatever length, and fill to tlen_bc with vpx_end
        - This is to match previous versions 'short' barcodes
        - I think we should be reporting the actual extracted barcode instead, but whatever!
4. fastq_parser.py is fundamentally different now, incl. new match_extract_dual_targ()

v2.40
1. Much time spent on making pyright happy (explicit type declarations, etc) - no functional changes

v2.41
1. Made output compatible with R-based PCRU concatination app (samp_date format problem)
2. Added 239X/INT, NIRM, OptM, and dGY stocks to config
3. Strip invisible whitespace from runinfo
4. Error check for known indexes - suggest runinfo format issue
5. Stress-test vs. large NIRM run

v2.42
1. Add columns to output runinfo:
    - n_reads, n_bc, n_reads_ac, n_bc_named, count-to-input_ratio (see io.py/write_output())
2. Modify output matrix:
    - Add sequencing input row to header
    - move bc_name and bc_seq to left

v2.43
1. added samp_date to matrix output
2. modified i/o error message in response to a single-vs-dual index issue w/ identical sample names.
3. modified to handle low input samples such that "empty" indexes with no barcodes don't break pipe. 

v2.50
1. Added modifiable quality params to adjust minimum req'd phred scores
    - mean_qual, mask_qual, and mismatches, passed through argparse or GUI
2. Changed date formats
3. Sort matrix output by date
4. New default run includes fuzzy matching and lowQ base masking
5. Changed argument handling across accessory scripts, all were reversioned
6. Substantially modified GUI
7. (Non-python) ./_misc_/MacroTemplate.xlsm included to recalculate proportions post-run.

v2.51
1. Moved unique/above cutoff to rowgroup 0 for excel sample sheets  
2. Recalculate rowgroup 0 props on ouput sample sheets to total 1

v2.52
1. Updated config to run all "M" stocks with dual-index
2. Changed fastq streaming behavior to report count rather than progress bar - saves a lot of time

v2.53
1. In cases where no named barcodes above cutoff are found, instead use nominal 1/input (min theoretical named) for uniques
2. Conditionally add per-index 'notes' column to runinfo
    - "ambig_children collapsed"
    - "no named above_cutoff"
3. Drop non-collapsed ambig uniques from results

Proposed:
Check for index compatability before counting?

Necessary:
Deal with nonspecific inputs! 0/Max fail.
"""
