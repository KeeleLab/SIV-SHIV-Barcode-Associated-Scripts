# io.py

"""
Name:       io.py
Author:     CAG
Version:    1.5
Date:       2025/6/04
"""

# %% Imports

import os
import csv
import gzip
import pandas as pd
from tqdm import tqdm
from typing import Generator
from utils.countdata_parser import parse_countdata
from utils.fastq_parser import parse_fastq_rec
from datetime import datetime

# %% Functions


def read_csv_to_dict(file_path: str):
    """
    Parse two-col csv into dict {'col1':'col1'}
    (Used for primer ref data, k=seq v=name)
    """
    result = {}
    with open(file_path, "r") as csvfile:
        reader = csv.reader(csvfile)
        # Skip the header row
        next(reader)
        for row in reader:
            col1, col2 = row
            result[col1] = col2
    return result


def read_fasta_to_dict(file_path: str):
    """
    Parse a standard fasta into a dict
    (Used for barcode ref data, k=name v=seq)
    """
    sequences = {}
    current_seq = []
    current_header = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                    current_seq = []
                current_header = line[1:]
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences


def format_refdata(runinfo_path: str, primer_path: str, barcode_path: str):
    """
    Generate reference data
    """
    print("Formatting reference data... \n")

    """     
    # Declare global vars to hold reference data
    global runinfo     # For runinfo.xlsx as dataframe
    global p5_refdict  # For dictionary of req'd p5 primer indexes as {name:seq}, p7 is optionally built in main()
    global bc_refdict  # For dictionary of barcodes as {name:seq} 
    """

    # Get runinfo
    runinfo = pd.read_excel(runinfo_path)
        # Strip whitespace from column names
    runinfo.columns = runinfo.columns.str.strip()
        # Strip whitespace from string data
    runinfo = runinfo.apply(lambda col: col.apply(lambda x: x.strip() if isinstance(x, str) else x))
        # Convert 'Date' column to datetime if not already, then extract just the date
    runinfo['Date'] = pd.to_datetime(runinfo['Date']).dt.date

    # Generate primer index reference as {name:seq}
    p5_refdict = read_csv_to_dict(primer_path)
    
    # Generate barcode reference as {name:seq}
    bc_refdict = read_fasta_to_dict(barcode_path)

    return (runinfo, p5_refdict, bc_refdict)


def stream_fastq(file_path: str) -> Generator[str, None, None]:
    """
    Iteratively stream fq or fq.gz entries into memory as 4-line string.
    Shows an indeterminate tqdm progress bar (no pre-counting).
    """
    open_func = gzip.open if file_path.endswith(".gz") else open

    with open_func(file_path, "rt") as fastq_file, tqdm(
        desc="Counting index/barcode pairs", unit="rec", ascii="-="
    ) as pbar:
        record = []
        for line in fastq_file:
            record.append(line.strip())
            if len(record) == 4:
                yield "\n".join(record) + "\n"
                record = []
                pbar.update(1)


def process_fastq_stream(
    sample_path: str, 
    settings: dict, 
    ):
    """
    Use stream_fastq() and parse_fastq_rec class to collect {idx:{bc:count}}
    """
    # Initialize dict to hold results
    res_pfq = {}

    for record in stream_fastq(sample_path):
        # split record to individual lines [header, seq, +, qual]
        fq = record.splitlines()

        # Parse record to return results.
        res = parse_fastq_rec(
            seq=fq[1],  # sequence
            qual=fq[3],  # qstring
            settings=settings,  # run settings
        )

        # Ignore read if barcode is missing or low mean q
        if not hasattr(res, "bc_seq") or res.bc_qual < settings["mean_qual"]:
            continue

        # Process p5 - skip if missing or low mean q
        if not (hasattr(res, "p5_seq") and hasattr(res, "p5_qual")):
            continue

        p5_seq = res.p5_seq if res.p5_qual >= settings["mean_qual"] else None
        if p5_seq is None:
            continue

        # Optionally process p7, fill instead of skip:
        # Ns where q<mask_qual, or all Ns if mean q<min_qual, X if missing
        p7_seq = ""

        if settings["dualindex"]:
            if hasattr(res, "p7_seq") and hasattr(res, "p7_qual"):
                p7_seq = res.p7_seq if res.p7_qual >= settings["mean_qual"] else "N" * settings["tlen_p7"]
            else:
                p7_seq = "X" * settings["tlen_p7"]

        # Set idx_seq as p5::p7, p7 is empty unless dualindex
        idx_seq = p5_seq + p7_seq

        # Add idx subdict to count dict if not present.
        if idx_seq not in res_pfq:
            res_pfq[idx_seq] = {}

        # Add bc to idx subdict if not present.
        if res.bc_seq not in res_pfq[idx_seq]:
            res_pfq[idx_seq][res.bc_seq] = 0

        # Increment counter for {idx:{bc:count}}.
        res_pfq[idx_seq][res.bc_seq] += 1

    return res_pfq


def write_output(
    runinfo_df: pd.DataFrame,
    settings_dict: dict,
    args_dict: dict,
    res_pcd: parse_countdata,
    out_path: str,
    legacy_format=False,
    filt_mat_ac=False,
):
    """
    Write the long data as <runname>_concat.csv, and a legacy excel file as <runname>_Analysis.csv
    Formatting is done via parse_countdata instance methods, see countdata_parser.py
    """
    
    def modify_runinfo(runinfo: pd.DataFrame, res_pcd: parse_countdata, missing: list) -> pd.DataFrame:
        """
        Modify runinfo DataFrame with additional analysis metadata
        
        Parameters:
        -----------
        runinfo : pd.DataFrame
            Original runinfo DataFrame
        res_pcd : parse_countdata
            Parsed count data object
        
        Returns:
        --------
        pd.DataFrame
            Modified runinfo DataFrame with additional columns
        """
        # Create a copy to avoid modifying the original DataFrame
        modified_runinfo = runinfo.copy()
        
        # =============== RUNINFO MODIFICATIONS ===============
        # Total Sequences Extracted per Primer (n bc detected w/ valid index)
        modified_runinfo['n_reads'] = modified_runinfo['full_idx'].map(
            {index: sum(barcode.values()) for index, barcode in res_pcd.idx_known.items()}
        )
        
        # Total distinct barcodes
        modified_runinfo['n_bc'] = modified_runinfo['full_idx'].map(
            {index: len(barcodes) for index, barcodes in res_pcd.idx_known.items()}
        )
        
        # Total named barcodes
        modified_runinfo['n_bc_named'] = modified_runinfo['full_idx'].map(
            res_pcd.df[~res_pcd.df["bc_name"].str.contains("Unique_", na=False)]
            .groupby("idx_name").size().to_dict()
        )
        
        # Total counts above cutoff
        modified_runinfo['n_reads_ac'] = modified_runinfo['full_idx'].map(
            res_pcd.df[(res_pcd.df["above_cutoff"] == True)]
            .groupby("idx_name")["bc_count"].sum().to_dict()
        )
        
        # Count to input ratio
        modified_runinfo['count-to-input_ratio'] = (
            modified_runinfo['n_reads_ac'] / modified_runinfo['Input TOTAL PER BARCODE']
        )

        # Add column if there are any notes
        if any(res_pcd.idx_notes.values()):
            modified_runinfo['notes'] = modified_runinfo['full_idx'].map(
                {index: '; '.join(notes) for index, notes in res_pcd.idx_notes.items() if notes}
            )

        return modified_runinfo

    # =============== END RUNINFO MODIFICATIONS ===============

    # Create the output directory if it doesn't exist
    os.makedirs(out_path, exist_ok=True)

    # Write long df
    res_pcd.df.to_csv(os.path.join(out_path, f"{res_pcd.rnam}_concat.csv"), index=False)

    # Obtain matrix and table outputs
    group_mats_dict = res_pcd.samp_group_matrices(filt_ac=filt_mat_ac)
    group_tabs_dict, _ = res_pcd.samp_group_sidelong_tables(legacy_format=legacy_format)

    # Get samp_groups from runinfo & output data - find/report missing (no barcodes!)
    expected_groups = runinfo_df["Animal"].dropna().unique()
    groups = res_pcd.df["samp_group"].unique()

    missing = set(expected_groups) - set(groups)

    if missing:
        print(f""" 
        No barcodes found for {missing}!
        """)

    # Check sample groups from runinfo against keys in res dicts
    if not set(group_mats_dict.keys()) == set(groups) or not set(group_tabs_dict.keys()) == set(groups):
        raise ValueError(f"""
        Check concat output for missing data. 
        Details:
        - group_mats_dict keys: {set(group_mats_dict.keys())}
        - group_tabs_dict keys: {set(group_tabs_dict.keys())}
        - expected groups: {set(groups)}

        Previously seen issues here:
        - ran 239M/M2 instead of 239M/M2_dual-index.
        """)

    # Modify runinfo with nested function
    runinfo_df = modify_runinfo(runinfo_df, res_pcd, missing)

    # Open excel writer for xlsx output
    with pd.ExcelWriter(os.path.join(out_path, f"{res_pcd.rnam}_Analysis.xlsx")) as writer:
        # Write the runinfo as sheet 1
        runinfo_df.to_excel(writer, sheet_name="runinfo", index=False)

        # For each samp_group...
        for samp_group in groups:
            print(f"writing {samp_group}...")

            # Write table
            group_tabs_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} samples", index=False, header=False
            )

            # Write matrix
            group_mats_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} matrix", index=True, header=True
            )

        # Add run metadata and arguments to settings_dict and create df
        settings_dict["analysis date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        settings_dict.update(args_dict)
        settings_df = pd.DataFrame(list(settings_dict.items()))

        # Add settings sheet to excel
        settings_df.to_excel(
            writer, sheet_name="Analysis settings", index=False, header=False
        )
# %%
