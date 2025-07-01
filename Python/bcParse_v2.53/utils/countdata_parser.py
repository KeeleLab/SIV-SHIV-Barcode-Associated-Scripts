# countdata_parser.py

"""
Name:       countdata_parser.py
Author:     CAG
Version:    1.5
Date:       2025/6/04
"""

# %% Imports

import pandas as pd
import numpy as np
import Levenshtein as lev
from tqdm import tqdm
from datetime import datetime
from collections import defaultdict

# %% Class definition

class parse_countdata(object):
    """
    Class to process the counts dictionary {idx_seq:{bc_seq:count}} into class object with attributes:
        - rawdict: {idx_seq:{bc_seq:count}}, raw counts
        - bc_refdict_kseq: {bc_seq:bc_name}, for mapping known barcodes
        - idx_expect_kseq: {idx_seq:idx_name}, for mapping for expected indexes
        - idx_known: {idx_name:{bc_seq:count}}, for named indexes
        - idx_other: {idx_seq:{bc_seq:count}}, for unknown indexes
        - bc_set: df with non-redundant barcode information
        - df: final processed dataframe

    output methods:
    parse_countdata.samp_group_matrices()
    parse_countdata.samp_group_sidelong_tables()
    """

    def __init__(
        self,
        countdict: dict,
        runinfo: pd.DataFrame,
        bc_refdict: dict,
        p5_refdict: dict,
        p7_refdict: dict | None,
    ):
        # Init `rawdict` {idx_seq:{bc_seq:count}}
        self.rawdict = countdict

        # Init `bc_refdict_kseq` {bc_seq:bc_name}
        self.bc_refdict_kseq = {value: key for key, value in bc_refdict.items()}

        # Init `idx_expect_kseq` {idx_seq:idx_name} for expected indexes
        # This will either create p5 only or p5::p7 version
        if p7_refdict is None:
            self.single_index_ref(p5_refdict, runinfo)
        else:
            self.dual_index_ref(p5_refdict, p7_refdict, runinfo)

        # Initialize `idx_known` {idx_name:[]}, empty dict with expected idx_name keys
        self.idx_known = {k: {} for k in self.idx_expect_kseq.values()}

        # Initialize `idx_other` {unknown_idx_seq:{bc_seq:count}}
        self.idx_other = {}

        # Populate `idx_known` and `idx_other` from `rawdict`
        self.combine_known_idx()

        # Generate `bc_set` and `df` from `idx_known`
        self.count_dict_to_dfs()

        # Add run meta to obj and merge runinfo to `df`
        self.add_runinfo(runinfo)

        # Per-index calculations and QC flags
        self.idx_notes = defaultdict(list)
        self.qc_ops()

    #############################
    # init methods
    #############################

    def single_index_ref(self, p5_refdict, runinfo):
        """
        Generate an index reference dictionary for p5-only runs, e/g {'VPX.P5.N':'ACGTACGT'}
        """
        self.idx_expect_kseq = {
            v: k for k, v in p5_refdict.items() if k in runinfo["Barcodes"].values
        }

    def dual_index_ref(self, p5_refdict, p7_refdict, runinfo):
        """
        Generate an index reference dictionary for dual index runs.
        This joins P5::P7 as a single sequence, e/g {'VPX.P5.N_VPX.P7':'ACGTACGTACGTACGT'}
        """
        # Concatenate p7 to p5
        self.idx_expect_kseq = {}

        # Pull only complete rows with these
        cleaned_runinfo = runinfo.dropna(subset=["Barcodes", "(F Barcode)"])

        # Iterate through each row in runinfo to build idx_expect_kseq
        for _, row in cleaned_runinfo.iterrows():
            idx_name = row["Barcodes"]
            p7_name = row["(F Barcode)"]

            # Check if index exists in p5_refdict
            if idx_name not in p5_refdict.keys():
                raise ValueError(f"'{idx_name}' from runinfo not in p5_refdict.")

            idx_seq = p5_refdict[idx_name]  # Get the corresponding index sequence

            # Check if index exists in p7_refdict
            if p7_name not in p7_refdict.keys():
                raise ValueError(f"'{p7_name}' from runinfo not in p7_refdict.")

            p7_seq = p7_refdict[p7_name]  # Get the corresponding p7 sequence

            # Concatenate sequences and create combined name
            full_seq = idx_seq + p7_seq
            full_name = f"{idx_name}_{p7_name}"

            # Add to the dictionary
            self.idx_expect_kseq[full_seq] = full_name

    def combine_known_idx(self):
        """
        Where raw idx seqs have exactly 1 ldist match to expected,
        combine per barcode dicts in idx_known dict {idx_name:{bc_seq:counts}}.
        No match or multimatch send to idx_other.
        """
        for k, _ in self.rawdict.items():
            # ldist each found idx against expected ref_idx seqs and get matches maxdist=1
            matches = self.ld_match_k_ret_v(k, self.idx_expect_kseq, 1)
            # if exactly one match, combine dict with idx_known using that match as the key
            if len(matches) == 1:
                self.idx_known[matches[0]] = {
                    key: self.idx_known[matches[0]].get(key, 0)
                    + self.rawdict[k].get(key, 0)
                    for key in set(self.idx_known[matches[0]]) | set(self.rawdict[k])
                }
            # else send to idx_other
            else:
                self.idx_other[k] = self.rawdict[k]
        
        # Check for known indexes
        if not any(self.idx_known.values()):
            raise ValueError(f"""
            Error! No expected primers found.
            Scroll up to find "primer_path" or "p7_path" for file location. 
            Check that runinfo columns match name format in primer reference fasta.
            """)

    def count_dict_to_dfs(self):
        """
        Convert dict = {'idx':{'bc':count}} to 3-col df to obtain a per-run set of named barcodes.
            - Uniques are named consistentily across all indexes in rank order.
            - Generates nonredundant df self.bc_set and complete self.named_df with counts.
        """
        # Generate the raw df
        rawdf = pd.DataFrame.from_records(
            [
                (idx_key, bc_key, bc_value)
                for idx_key, idx_dict in self.idx_known.items()
                for bc_key, bc_value in idx_dict.items()
            ],
            columns=["idx_name", "bc_seq", "bc_count"],
        )

        # Dealing with uniques:
        # Get df with total counts per barcode sequence
        df = (
            rawdf.groupby(["bc_seq"])["bc_count"]
            .sum()
            .reset_index()
            .sort_values(by="bc_count", ascending=False)
        )

        # map barcodes names from refdict
        df["bc_name"] = df["bc_seq"].map(self.bc_refdict_kseq)

        # swap unmatched NaN for 'Unique_'
        df["bc_name"] = df["bc_name"].replace(np.nan, "Unique_")

        # Identify duplicate 'bc_name' values, which should only include 'Unique_'
        mask = df["bc_name"].duplicated(keep=False)

        # Add a counter to duplicate 'Unique_' values
        df.loc[mask, "bc_name"] = self.add_unique_counter(df.loc[mask])

        # Sort by count and select seq/name columns
        df = df.sort_values(by="bc_count", ascending=False)
        df = df[["bc_seq", "bc_name"]]

        # Save a complete set of named barcode seqs for the run
        self.bc_set = df.reset_index(drop=True)

        # Merge df to rawdf to reassign names
        df = pd.merge(rawdf, self.bc_set, how="left", on="bc_seq")
        df = df.sort_values(["idx_name", "bc_count"], ascending=[True, False])

        # If barcode name is found twice, it's in multiple indexes.
        df["x_idx"] = df.duplicated(subset=["bc_name"], keep=False)

        # Set self.named_df
        self.df = df.reset_index(drop=True)

    def add_runinfo(self, runinfo):
        """
        Clean and rename runinfo columns, merge to named_df
        """
        # Get run metadata as object attrs
        self.rnum = runinfo["Run Number"].iloc[0]

        self.rnam = runinfo.loc[
            runinfo.index[runinfo["Run Number"] == "Run Name"].tolist()[0] + 1,
            "Run Number",
        ]

        self.rdat = runinfo.loc[
            runinfo.index[runinfo["Run Number"] == "Date"].tolist()[0] + 1, "Run Number"
        ]

        # Subset and rename per-idx metadata columns from runinfo
        self.idx_meta = runinfo.dropna(subset=["Barcodes", "(F Barcode)"])[
            ["Animal", "Sample", "full_idx", "Date", "Input TOTAL PER BARCODE"]
        ].rename(
            columns={
                "Animal": "samp_group",
                "Sample": "samp_name",
                "full_idx": "idx_name",
                "Date": "samp_date",
                "Input TOTAL PER BARCODE": "input",
            }
        )

        # Merge runinfo to named_df and raise error if any rows are missed
        try:
            merged_df = pd.merge(self.df, self.idx_meta, on="idx_name", how="left")

            # Check if any rows in self.named_df didn't find a match
            if merged_df["samp_group"].isnull().any():
                raise ValueError("rows in named_df did not find idx match in runinfo")

            # Update self.named_df with the merged result
            self.df = merged_df

        except ValueError as e:
            # Re-raise the ValueError with additional context if needed
            raise ValueError(f"Merge operation failed: {str(e)}") from e

    def qc_ops(self):
        """
        Perform QC operations:

        - Global
            - flag x_group barcodes (putative idx_hops)
            - flag short barcodes
            - final column arrangement

        - Index specific
            - calculate proportions
            - flag above-cutoff
            - ldist check for putative parents

        This results in a complete, 'final' long-format dataframe
        experimental: self.idx_notes are built for io.py runinfo updates

        See class & static methods below.
        """
        # flag cross-group barcodes
        self.df["x_group"] = (
            self.df.groupby("bc_name")["samp_group"].transform("nunique") > 1
        )

        # Split df by idx_name to {'idx_name':pd.DataFrame} dict to process individually
        df_dict = {name: group for name, group in self.df.groupby("idx_name")}

        # Per-index processing
        for idx, df in df_dict.items():
            
            # proportions
            df = self.get_prop(df)

            # above cutoff
            df, note = self.flag_above_cutoff(df)
            note is not None and self.idx_notes[idx].append(note)

            # ldist/putative parents
            df = self.flag_putative_parents(idx, df)

            # collapse ambiguous children
            df, note = self.collapse_ambig_children(df)
            note is not None and self.idx_notes[idx].append(note)

            df_dict[idx] = df

        # Concatenate the modified DataFrames
        self.df = pd.concat(df_dict.values(), axis=0)

        # Flag short barcodes
        self.df = self.flag_short(self.df)

        # Arrange columns and rows
        self.df = self.df.loc[
            :,
            [
                "idx_name",  # Index name - P5 or P5_P7
                "samp_group",  # Sample Group - animal, donor, etc
                "samp_name",  # Sample Name - specific sample
                "samp_date",  # From runinfo
                "input",  # ""
                "bc_name",  # Barcode name
                "bc_count",  # Barcode miseq counts
                "proportion",  # Per-index proportion
                "above_cutoff",  # Proportion > 1/input (named) or 100 * min named prop (unique)
                "x_idx",  # Cross-index barcode?
                "x_group",  # Cross-group barcode?
                "short_bc",  # Short barcode?
                "putative_parent",  # ldist hit to named barcode in this index?
                "bc_seq",  # Barcode sequence
            ],
        ].sort_values(["idx_name", "bc_count"], ascending=[True, False])

        # Add run metadata as columns
        add_cols = {"run_num": self.rnum, "run_name": self.rnam, "run_date": self.rdat}

        self.df = pd.concat(
            [pd.DataFrame(add_cols, index=self.df.index), self.df], axis=1
        )

        # Reset index
        self.df.reset_index(drop=True, inplace=True)

    def flag_putative_parents(self, idx_name: str, df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()

        Necessarily a class method to allow internal calls to static methods - it wasn't working as static, anyway.

        Creates a column in self.df that reports highest count barcode with ldist = 1 to current row.
        Runs by nesting static methods get_ldist_matrix, get_countdif_matrix within get_filtered_hits.
        """
        print(f"For {idx_name}: ")

        # Subset unique and known
        df_unique = df[df["bc_name"].str.startswith("Unique_")]
        df_known = df[~df["bc_name"].str.startswith("Unique_")]

        # Get parents for known vs. known
        known_parents = self.get_filtered_hits(
            self.get_ldist_matrix(  # Get known ldist matrix
                df_known["bc_seq"].values,
                df_known["bc_seq"].values,
                tqdm_desc="  Calculating known ldists",
            ),
            self.get_countdif_matrix(  # Get known countdiff matrix
                df_known["bc_count"].values,
                df_known["bc_count"].values,
                tqdm_desc="  Resolving known parents",
            ),
            # Pass known barcode names in twice
            df_known["bc_name"].values,
            df_known["bc_name"].values,
        )

        # Get parents for unique vs. known
        unique_parents = self.get_filtered_hits(
            self.get_ldist_matrix(  # Get unique vs known ldist matrix
                df_unique["bc_seq"].values,
                df_known["bc_seq"].values,
                tqdm_desc="  Calculating unique ldists",
            ),
            self.get_countdif_matrix(  # Get unique vs known countdiff matrix
                df_unique["bc_count"].values,
                df_known["bc_count"].values,
                tqdm_desc="  Resolving unique parents",
            ),
            # Pass unique and known barcode names
            df_unique["bc_name"].values,
            df_known["bc_name"].values,
        )

        all_parents = known_parents | unique_parents
        df["putative_parent"] = df["bc_name"].map(all_parents)

        print("\n")

        return df

    @staticmethod
    def ld_match_k_ret_v(word: str, dictionary: dict, max_distance: int):
        """
        Get LD of word vs. dict keys, return values for keys below max as list
        """
        matches = []
        for key in dictionary:
            distance = lev.distance(word, key, score_cutoff=2)
            if distance <= max_distance:
                matches.append(dictionary[key])
        return matches

    @staticmethod
    def add_unique_counter(df: pd.DataFrame):
        """
        'Uniqify' barcode names by adding number observed, in order of count.
        This will only affect replicate barcode names (in this case, looking at totals per run, 'Unique' is only affected.)
        """
        return df["bc_name"] + df["bc_name"].groupby(df["bc_name"]).cumcount().add(
            1
        ).astype(str)

    @staticmethod
    def get_prop(df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()
        """
        # calculate proportion
        df["proportion"] = df["bc_count"] / df["bc_count"].sum()

        return df

    @staticmethod
    def flag_above_cutoff(df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()

        Assign above_cutoff as True for rows where:
            - bc_name does not contain 'Unique_' and proportion is greater than 1/input
            - bc_name contains 'Unique' and bc_count is greater than minimum named proportion * 100
        """
        # Initialize above_cutoff col
        df["above_cutoff"] = False

        # Set named cutoff
        named_cutoff = 1 / df["input"].iloc[0]

        # Flag above_cutoff for non-unique barcodes
        df.loc[
            ~df["bc_name"].str.contains("Unique_")
            & (df["proportion"] > named_cutoff),
            "above_cutoff",
        ] = True

        # Get minimum named above_cutoff proportion
        min_named_p = df.loc[df["above_cutoff"], "proportion"].min()

        # Substitute nominal value if min_named_p is missing/nan
        if np.isnan(min_named_p):
            unique_cutoff = 1 / df["input"].iloc[0]
            note = f"No named above_cutoff: used prop > {unique_cutoff:.3g}"
        # Otherwise use 100X minimum named proportion
        else: 
            unique_cutoff = 100 * min_named_p
            note = None

        # Flag above_cutoff for unique barcodes
        df.loc[
            df["bc_name"].str.contains("Unique_")
            & (df["proportion"] > unique_cutoff),
            "above_cutoff",
        ] = True

        return df, note

    @staticmethod
    def collapse_ambig_children(df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()

        Collapse children w/ single N to valid parent rows, and drop multi-N
        """
        # Skip if no Ns detected (as with unmasked runs, etc)
        if not df["bc_seq"].str.contains("N").any():
            note = None
            return df, note

        note = "ambig_children collapsed"

        # Drop rows where 'bc_seq' contains more than one 'N'
        df = df[~(df["bc_seq"].str.count("N") > 1)]

        # Identify ambiguous children
        ambig_children = df[
            (df["putative_parent"].notna())  # has named putative parent
            & (
                df["bc_seq"].str.contains("N")
            )  # contains a single N (multi was dropped)
            & (~df["above_cutoff"])  # is below cutoff
            & (
                df["putative_parent"].isin(df[df["above_cutoff"]]["bc_name"])
            )  # parent is above cutoff
        ]

        # Drop these (and remaining single-Ns) from df
        #df = df.drop(ambig_children.index)
        df = df[~(df["bc_seq"].str.count("N") > 0)]

        # Collapse ambig child counts by parent
        add_counts = ambig_children.groupby("putative_parent")["bc_count"].sum()

        # Add counts to parent rows in df
        df.loc[df["bc_name"].isin(add_counts.index), "bc_count"] += (
            df["bc_name"].map(add_counts).fillna(0)
        )

        # Recalculate the proportion column
        tot = df["bc_count"].sum()
        df["proportion"] = df["bc_count"] / tot if tot > 0 else 0

        return df, note

    @staticmethod
    def get_ldist_matrix(seqs_A: list[str], seqs_B: list[str], tqdm_desc: str):
        """
        Generic method to return a levenshtein distance matrix of barcode sequence sets A and B
        """
        # Initialize and populate levenshtein distance matrix
        distmat = np.zeros((len(seqs_A), len(seqs_B)), dtype=int)

        for i, A in enumerate(tqdm(seqs_A, desc=tqdm_desc, ascii="-=")):
            for j, B in enumerate(seqs_B):
                distmat[i, j] = lev.distance(A, B)

        return distmat

    @staticmethod
    def get_countdif_matrix(counts_A: list[int], counts_B: list[int], tqdm_desc: str):
        """
        Generic method to return a matrix of barcode count differences for sets A and B
        """
        # Initialize and populate count-difference matrix
        diffmat = np.zeros((len(counts_A), len(counts_B)), dtype=int)

        for i, A in enumerate(tqdm(counts_A, desc=tqdm_desc, ascii="-=")):
            for j, B in enumerate(counts_B):
                diffmat[i, j] = B - A

        return diffmat

    @staticmethod
    def get_filtered_hits(
        distmat: np.ndarray, diffmat: np.ndarray, names_A: list[str], names_B: list[str]
    ):
        """
        Return dictionary of {child, parent} pairs where:
            - max ldist is 1
            - putative parent has the greatest count difference to child
        In ties, it will return the first value.
        """
        # Set a mask for desired distance threshold
        mask = distmat == 1

        # Return a dict of {A, B} values where known had the highest count
        results = {}
        for i, A in enumerate(names_A):
            # Get matrix indexes where distance == threshold
            matches = np.where(mask[i])[0]
            # If they exist...
            if len(matches) > 0:
                # Get the index of the named match with the highest count difference
                j = matches[np.argmax(diffmat[i, matches])]
                # If the difference is positive (B > A, thus putative parent):
                if diffmat[i, j] > 0:
                    # Use index to query the name and update dict
                    results[A] = names_B[j]

        return results

    @staticmethod
    def flag_short(df: pd.DataFrame):
        """
        Simple method to flag short barcodes
        Currently catches 10 or more bases of VPX
        This could be improved by dynamically getting 'fill_seq' end
        OR
        passing short status forward with fastq_parser
        """
        # Initialize column
        df["short_bc"] = False
        # Flag true if motif detected
        df["short_bc"] = df["bc_seq"].str.contains("ACTAGCATAA")

        return df

    ###################
    # external methods
    ###################

    def samp_group_matrices(self, filt_ac=False):
        """
        Return a per-samp_group dict of barcode proportion matrices,
        optionally filtered for above_cutoff
        """
        df = self.df.copy()

        # Get distinct {samp_group: [idx_name(s)]}
        groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

        # Initialize dict to hold per-samp_group results
        df_mats = {}

        for samp_group, idx_names in groups.items():
            
            # Optional cutoff filter
            if filt_ac:
                # filter above cutoff rows
                df_filt = df[(df["idx_name"].isin(idx_names) & df["above_cutoff"])].copy()
                # recalculate proportions based on above_cutoff data
                df_filt["proportion"] = (
                    df_filt.groupby("idx_name")["proportion"]
                    .transform(lambda x: x / x.sum())
                )
            else:
                df_filt = df[(df["idx_name"].isin(idx_names))]

            # Generate bc proportion matrix
            df_mat = (
                df_filt.pivot(
                    index=["bc_name", "bc_seq"],
                    columns=["samp_name", "samp_date", "idx_name", "input"],
                    values="proportion",
                )
                .fillna(0)
                .sort_index(axis=1, level=["samp_name", "samp_date"])
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

        return df_mats

    def samp_group_sidelong_tables(self, legacy_format=False):
        """
        "Sidelong tables" are generated for each samp_group and held in a dict.

        This generates a viewer-oriented table containing adjascent long dataframes
        which indicate cutoffs via padding rows. It's not computer-friendly.

        The legacy format includes a weird header and the original colnames, and
        (I hate that) it's included to ensure continuity with legacy data.

        Really struggled in making this neat, elegant, or modular.
        """
        df = self.df.copy()

        # Get samp_group indexes
        groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

        # Split long df to dict w/ keys = idx_name
        split_dfs = {k: v for k, v in df.groupby("idx_name")}

        # Per-index data formatting
        for idx in split_dfs:
            # subset index's df
            df = split_dfs[idx].copy()

            #########
            # body df
            #########

            df_body = df.copy()[
                [
                    "bc_name",
                    "bc_seq",
                    "bc_count",
                    "proportion",
                    "putative_parent",
                    "x_group",
                    "short_bc",
                    "above_cutoff",
                ]
            ]

            # Conditions for grouping barcodes
            c1 = df_body["above_cutoff"]  # Above cutoff
            c2 = df_body["putative_parent"].notna()  # Has putative parent
            c3 = df_body["bc_name"].str.contains("Unique")  # Is Unique

            # Bitpack the flag columns - this is actually slick.
            df_body["flag"] = (
                c1.astype(int) * 4 + c2.astype(int) * 2 + c3.astype(int) * 1
            )

            """ 
            # Possible flag values:
            flags = {
                0: 'Below cutoff, No parent, Named',
                1: 'Below cutoff, No parent, Unique',
                2: 'Below cutoff, Has parent, Named',
                3: 'Below cutoff, Has parent, Unique',
                4: 'Above cutoff, No parent, Named',
                5: 'Above cutoff, No parent, Unique',
                6: 'Above cutoff, Has parent, Named',
                7: 'Above cutoff, Has parent, Unique'
            } 

            # With this format, the groups we want are:
            #   0: 4,5
            #   1: 6,7
            #   2: 0,1,2,3 (the rest)
            """

            # map the grouping variables
            df_body["rowgroup"] = df_body["flag"].map({4: 0, 5: 0, 6: 1, 7: 1}).fillna(2)

            # sort by ascending group and descending proportion
            df_body = df_body.sort_values(
                by=["rowgroup", "proportion"], ascending=[True, False]
            ).reset_index(drop=True)

            # I'm pained by this: 
            # Recalculate group 0's proportions to equal 1...
            ac_tot_prop = df_body.loc[df_body['rowgroup'] == 0, 'proportion'].sum()

            if ac_tot_prop != 0:
                df_body.loc[df_body['rowgroup'] == 0, 'proportion'] /= ac_tot_prop
            else:
                # Handle edge case if sum is zero
                df_body.loc[df_body['rowgroup'] == 0, 'proportion'] = 0

            # find indices where the rowgroup changes
            changes = df_body.index[
                df_body["rowgroup"] != df_body["rowgroup"].shift()
            ].tolist()
            changes.pop(
                0
            )  # artifact of using shift() causes 'change' at first row - pop this.

            # insert empty rows
            pad = pd.DataFrame(
                [[np.nan] * len(df_body.columns)], columns=df_body.columns
            )

            for index in sorted(changes, reverse=True):
                df_body = pd.concat(
                    [df_body.iloc[:index], pad, df_body.iloc[index:]], ignore_index=True
                )

            # drop columns
            df_body = df_body.drop(["flag", "rowgroup", "above_cutoff"], axis=1)

            # Format-specific modification
            if legacy_format:
                # Cast these to string for subsequent map operation
                df_body[["x_group", "short_bc"]] = df_body[
                    ["x_group", "short_bc"]
                ].astype(str)

                # Get mask for spaced rows
                nan_mask = ~df_body["bc_name"].isna()

                # modify (string-cast) boolean columns
                df_body.loc[nan_mask, ["x_group", "short_bc"]] = df_body.loc[
                    nan_mask, ["x_group", "short_bc"]
                ].map(lambda x: "Yes" if x == "1.0" else "" if x == "0.0" else "")

                df_body.loc[~nan_mask, ["x_group", "short_bc"]] = df_body.loc[
                    ~nan_mask, ["x_group", "short_bc"]
                ].map(lambda x: "" if x == "nan" else x)

                # modify putative parent column
                df_body.loc[nan_mask, "putative_parent"] = (
                    df_body.loc[nan_mask, "putative_parent"]
                    .fillna("NA-NA")
                    .map(lambda x: f"{x}-1" if x != "NA-NA" else x)
                )

                # modify column names
                df_body = df_body.rename(
                    columns={
                        "bc_name": "Barcode",
                        "bc_seq": "Sequence",
                        "bc_count": "Counts",
                        "proportion": "Proportion",
                        "putative_parent": "Hamming Dist to Another Sequence",
                        "x_group": "Multi-group?",
                        "short_bc": "Short barcode?",
                    }
                )

            else:
                # modify boolean columns
                df_body[["x_group", "short_bc"]] = df_body[["x_group", "short_bc"]].map(
                    lambda x: True if x == 1 else False if x == 0 else ""
                )

            # Further padding steps for pd.concat to header:
            # - get column range body + 1
            new_cols = [f"{i}" for i in range(len(df_body.columns) + 1)]

            # - add original colnames as row
            df_body = pd.concat(
                [pd.DataFrame([df_body.columns], columns=df_body.columns), df_body]
            ).reset_index(drop=True)

            # - add padding columns
            for i in range(len(new_cols) - len(df_body.columns)):
                df_body[i] = np.nan

            # - Rename columns for concat
            df_body.columns = new_cols

            ###########
            # header df
            ###########

            # Initialize dict to hold final header values
            header_dict = {}

            # Check these columns for single value
            header_keys = [
                "run_num",  # run-specific
                "run_name",  # ""
                "run_date",  # ""
                "idx_name",  # Index-specific
                "samp_group",  # ""
                "samp_name",  # ""
                "samp_date",  # ""
                "input",  # ""
            ]

            for k in header_keys:
                v = df[k].unique()
                if len(v) != 1:
                    raise ValueError(f"{k} has multiple values: {v}")

                # If unique, update dict
                header_dict[k] = v[0]

            # Format-specific modification
            if legacy_format:

                # Modify date format
                header_dict["samp_date"] = header_dict["samp_date"].strftime('%Y-%m-%d')

                header_data = [
                    # Weird repeated sample metadata
                    [
                        header_dict["samp_group"],
                        header_dict["samp_name"],
                        header_dict["samp_date"],
                        header_dict["samp_group"],
                        header_dict["samp_name"],
                        header_dict["samp_date"],
                    ]
                    + [np.nan] * (len(new_cols) - 6),
                    # Actual "column names"
                    [
                        "Sequencer Input",
                        "Total # Sequences",
                        "# of Barcodes (Excluding Uniques)",
                        "Proportion Uniques",
                    ]
                    + [np.nan] * (len(new_cols) - 4),
                    # Values
                    [
                        header_dict["input"],
                        df["bc_count"].sum(),
                        df[~df["bc_name"].str.contains("Unique")].shape[0],
                        df[df["bc_name"].str.contains("Unique")]["bc_count"].sum()
                        / df["bc_count"].sum(),
                    ]
                    + [np.nan] * (len(new_cols) - 4),
                    # pad
                    [np.nan] * len(new_cols),
                ]

                df_head = pd.DataFrame(header_data, columns=new_cols)

            else:
                # Create the DataFrame with keys and values
                df_head = pd.DataFrame(list(header_dict.items()), columns=new_cols[:2])

                # pad width
                for col in new_cols[2:]:  # Start from the third column of new_cols
                    df_head[col] = np.nan

                # pad length
                pad = pd.DataFrame([[np.nan] * len(new_cols)], columns=new_cols)
                df_head = pd.concat([df_head, pad], ignore_index=True)

                # pad thai?

            # IT'S ALIVE - concat head and body
            df = pd.concat([df_head, df_body], ignore_index=True)

            # Update original dict with final df
            split_dfs[idx] = df

        # Human Centipede - join these sidelong
        sidelongs = {}

        for k, v in groups.items():
            # Extract DataFrames for each key
            dfs = [split_dfs.get(k2) for k2 in v if k2 in split_dfs]

            # Join horizontally and save to new dict
            sidelongs[k] = pd.concat(dfs, axis=1)

        return sidelongs, split_dfs