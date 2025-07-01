

#%% 

#sample_path = "../../_sampledata/Example/ExampleSetup_testing.fastq.gz"
#runinfo_path = "../../_sampledata/Example/ExampleSetupRunInfo.xlsx"
#out_path = "../../_sampledata/Example"
#stock = "239M"
#dualindex = False

sample_path = "../../_sampledata/P7newIndex2test/P7newIndex2-test-011325_S1_L001_R1_001.fastq.gz"
runinfo_path = "../../_sampledata/P7newIndex2test/P7new.Index2 test miseq Run Info Form.xlsx"
out_path = "../../_sampledata/P7newIndex2test"

#sample_path = "/Volumes/ACVP-RES-VEC-FS/Abbey/BaseSpace/Kulpa Runs/Kulpa LNMC/Kulpa-LNMC-040825_S1_L001_R1_001.fastq"
#runinfo_path = "/Volumes/ACVP-RES-VEC-FS/Abbey/BaseSpace/Kulpa Runs/Kulpa LNMC/Kulpa LNMC DNA Run 3 Miseq Info Sheet.xlsx"
#out_path = "/Users/goodmanca/Desktop/scratch/bcParse_troubleshooting/"

stock = "239M_dual-index"
#stock = "239M2"
mask = False
filt_mat_ac = True
legacy_format = True

mean_qual = 0
mismatches = 1
mask_qual = 20

#%%

# Load modules from bcParse.py!

# Get and print run settings - add new stocks and settings in config.py!
settings = get_settings(stock)
print("\n".join(f"{key}: {value}" for key, value in settings.items()) + "\n")

"""
2. Stream fastq to count barcodes reads per index

- Return {idx:{bc:count}}
- See utils/io.py and utils/fastq_parser.py
- mask returns N for base q<30
"""
res_pfq = process_fastq_stream(
    sample_path=sample_path, 
    settings=settings,
    mean_qual=mean_qual,
    mismatches=mismatches,
    mask=mask,
    mask_qual=mask_qual
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

#%%

""" 
4. Main output - see io.py & utils/countdata_parser.py. 
"""
# Collect args in dict to pass to output
args_dict = {
    "sample_path": sample_path,
    "runinfo_path": runinfo_path,
    "out_path": out_path,
    "stock": stock,
    "mask": mask,
    "filt_mat_ac": filt_mat_ac,
    "legacy_format": legacy_format,
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


#%%

a,b = res_pcd.samp_group_sidelong_tables(legacy_format=True)
c = res_pcd.samp_group_matrices(filt_ac=True)

df = res_pcd.df.copy()

# Get distinct {samp_group: [idx_name(s)]}
groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

# Initialize dict to hold per-samp_group results
df_mats = {}

filt_ac = True

for samp_group, idx_names in groups.items():
    # Optional cutoff filter
    if filt_ac:
        df_filt = df[(df["idx_name"].isin(idx_names) & df["above_cutoff"])]
    else:
        df_filt = df[(df["idx_name"].isin(idx_names))]

    df_mats[samp_group] = df_filt

    # Generate bc proportion matrix
    df_mat = (
        df_filt.pivot(
            index=["bc_name", "bc_seq"],
            columns=["samp_name", "idx_name", "input"],
            values="proportion",
        )
        .fillna(0)
        .reset_index()
    )

    # Get all columns except the first two
    cols_to_sort = df_mat.columns[2:]

    # Sort the DataFrame
    df_mat = df_mat.sort_values(
        by=list(cols_to_sort),
        ascending=[False] * len(cols_to_sort),
        na_position="last",
    ).reset_index(drop=True)

    # Move 'bc_seq' and 'bc_name' to the end
    #df_mat["bc_name"] = df_mat.pop("bc_name")
    #df_mat["bc_seq"] = df_mat.pop("bc_seq")

    # Add to dict
    df_mats[samp_group] = df_mat


#%%

""" 
Optional subanalyis for dualindex runs
"""
if settings['dualindex']:

    print(f"\n")
    print(f"Processing ambigous p7 reads...")
    print(f"\n")
    
    # Obtain a pfq dict of valid p5 but ambig p7
    res_pfq_ambig_p7 = {
        key[:-settings['tlen_p7']]: value 
        for key, value in res_pcd.idx_other.items()
        if set(key[:settings['tlen_p5']]).issubset('ACGT')
    }

    # Reset this column in runinfo since we've already written output
    runinfo['full_idx'] = runinfo['Barcodes']

    # Run pcd on subset
    res_pcd_ambig_p7 = parse_countdata(
        countdict=res_pfq_ambig_p7,
        runinfo=runinfo,
        bc_refdict=bc_refdict,
        p5_refdict=p5_refdict,
        p7_refdict=None
    )

    # Write df
    res_pcd_ambig_p7.df.to_csv(
        os.path.join(out_path, f"{res_pcd_ambig_p7.rnam}_concat_ambig_p7.csv"), index=False)

#%%

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

    # Add a new column to mark rows where 'Animal' is in 'missing'
    modified_runinfo['no_barcodes'] = (
        (~modified_runinfo['full_idx'].isnull()) &  # Ensure 'full_idx' is not NaN
        (modified_runinfo['Animal'].isin(missing))  # Check if 'Animal' is in 'missing'
    )
    
    return modified_runinfo


# %%

runinfocopy = runinfo.copy()

missing = {}

runinfo_mod = modify_runinfo(runinfocopy, res_pcd, missing)
# %%
