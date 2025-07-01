# config.py

"""
Name:       config.py
Author:     CAG
Version:    1.2
Date:       2025/2/25
"""

# %% Imports

import os

# %% User-defined settings

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the path to the reference data directory relative to the script
ref_data_dir = os.path.join(script_dir, "ref")

# Define stock settings
"""
General read structures - 

VPX runs; 
239M/M2/dual-index, NIRM, OptM, dGYM: 
    < 4bp tag :: (`tlen_p5`)bp idx :: `ref_p5`+ :: revcomp `ref_bc` :: revcomp (`tlen_bc`)bp barcode  >

    In 239M/M2, the biological data is reverse-complemented, but the index isn't, so the algotithm flow is:
    1. Identify `ref_p5` on raw sense strand
    2. Pull `tlen_p5` bases `upstream` of that
    3. `rdir` is rev, so reverse complement the sequence
    4. Identify `ref_bc` on revcomped new sense strand
    5. Pull `tlen_bc` bases `upstream` from the revcomped new sense strand

    In general, the index will always be apparent on the 5` end of the raw read.
    Revcomp occurs after identifying the index, when necessary for pulling region of interest.

    When running dualindex - modified read structure is as follows:
    < 239M/M2 read as above > :: < revcomp `ref_p7` :: revcomp (`tlen_p7`)bp p7 >

    NOTE that to match short barcodes from 239M/M2 p5-only characterization, provide 
    a 'vpx_end' sequence which automatically pads barcodes shorter than 34bp as  

    < 34 - length(bc_seq) bases 5' vpx_end 3' :: 5' bc_seq 3' > 

    ... see fastq_parser/match_extract_dual_targ(). 
    If vpx_end is not provided, the actual short sequence will be returned.

239X/INT:
    < 4bp tag :: (`tlen_p5`)bp idx :: `ref_p5` = `ref_bc` :: (`tlen_bc`)bp barcode ... >

    More straightforward. Read is fwd-oriented, and the p5 index and barcode flank a common SIV region.
    Like 239M/M2, the only difference is barcode reference data.
"""

settings_opts = {
    "239M_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "239M_reference.fasta"),
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },
    "239M2_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "239M2_reference.fasta"),
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },
    "OptM_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "OptM_reference.fasta"),
        "dualindex": False,
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },

    "NIRM_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "NIRM_reference.fasta"),
        "dualindex": False,
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },

    "dGYM_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "dGYM_reference.fasta"),
        "dualindex": False,
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },
    "V67M_dual-index": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "V67M_reference.fasta"),
        "dualindex": False,
        # dualindex settings
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": os.path.join(ref_data_dir, "P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },
    "239X/INT": {
        "rdir": "fwd",
        "ref_p5": "TAAAAATTTTCGGGTCTATTAC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "TAAAAATTTTCGGGTCTATTAC",
        "tdir_bc": "downstream",
        "tlen_bc": 45,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "239X_reference.fasta"),
        "dualindex": False,
    },
    "239M": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "239M_reference.fasta"),
        "dualindex": False,
    },
    "239M2": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "239M2_reference.fasta"),
        "dualindex": False,
    },
    "OptM": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "OptM_reference.fasta"),
        "dualindex": False,
    },
    "NIRM": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "NIRM_reference.fasta"),
        "dualindex": False,
    },
    "dGYM": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "dGYM_reference.fasta"),
        "dualindex": False,
    },
    "V67M": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": os.path.join(ref_data_dir, "P5_primers.csv"),
        "barcode_path": os.path.join(ref_data_dir, "V67M_reference.fasta"),
        "dualindex": False,
    },
}

# List of available stocks
stocks = list(settings_opts.keys()) + ["Manual CLI"]

# %% Functions


# Define stock dict
def get_settings(stock: str):
    if stock in settings_opts:
        settings_dict = settings_opts[stock].copy()
    else:
        settings_dict = {
            "rdir": input("Enter read direction: "),
            "ref_p5": input("Enter p5 reference sequence: "),
            "tdir_p5": input("Enter p5 target direction: "),
            "tlen_p5": int(input("Enter p5 target length: ")),
            "ref_bc": input("Enter barcode reference sequence: "),
            "tdir_bc": input("Enter barcode target direction: "),
            "tlen_bc": int(input("Enter barcode target length: ")),
            "primer_path": input("Path/to/sequencing_primers.csv"),
            "barcode_path": input("Path/to/barcode_reference.fasta"),
            "dualindex": input("True/False to run dual index mode"),
            "ref_p7": input("Enter p7 reference sequence: "),
            "tdir_p7": input("Enter p7 target direction: "),
            "tlen_p7": int(input("Enter p7 target length: ")),
            "p7_path": input("Path/to/p7_primers.csv"),
            "fill_seq": input("Pad sequence for dual-index 3' repair - ask Charlie"),
        }

    return settings_dict
