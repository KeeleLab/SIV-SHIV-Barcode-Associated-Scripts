# fastq_parser.py

"""
Name:       fastq_parser.py
Author:     CAG
Version:    1.3
Date:       2025/4/14
"""

# %% Imports

#import re
import regex
import numpy as np

# %% Class definition


class parse_fastq_rec(object):
    """
    Class to handle each fastq record as an instance while it streams.
    Return per-read idx (named if expected), bc, and assoc mean qualities.
    """

    def __init__(
        self, 
        seq: str, 
        qual: str, 
        settings: dict, 
    ):
        # Initialize sequence and quality
        self.seq = seq
        self.qual = self.Q33conv(qual)

        # (Opt) mask lowQ bases
        if settings["mask"]:
            self.seq = "".join(
                "N" if q < settings["mask_qual"] else b for b, q in zip(self.seq, self.qual)
            )

        # Process the fastq record according to settings
        self.process_sequence(settings)

        # Get read mean qual
        self.read_q = np.mean(self.qual)

    def process_sequence(self, settings):
        """
        Extract regions of interest and assoc. quality strings and update class attrs.
        Final attrs:
        - seq, qual
        - p5_seq, p5_qual
        - bc_seq, bc_qual
        - (Opt) p7_seq, p7_qual
        """
        # Extract p5 index
        p5 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            ref=settings["ref_p5"],
            targ_dir=settings["tdir_p5"],
            targ_len=settings["tlen_p5"],
            mismatches=settings["mismatches"],
        )

        # If p5 found, set attrs and continue
        if p5[0]:
            self.p5_seq, self.p5_qual = p5[0], np.mean(p5[1])

            # Reverse complement if called
            if settings["rdir"] == "rev":
                self.seq = self.revcomp(self.seq)
                self.qual = self.qual[::-1]

            # Process as single or dual index
            if not settings["dualindex"]:
                self.process_single_index(settings)
            else:
                self.process_dual_index(settings)

    def process_single_index(self, settings):
        """
        Obtain barcode based on target and length
        """
        # Extract barcode for single index
        bc = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            ref=settings["ref_bc"],
            targ_dir=settings["tdir_bc"],
            targ_len=settings["tlen_bc"],
            mismatches=settings["mismatches"],
        )

        if bc[0]:
            self.bc_seq, self.bc_qual = bc[0], np.mean(bc[1])

    def process_dual_index(self, settings):
        """
        Obtain p7, and if present obtain barcode based on two targets to 
        facilitate optional repair of short barcode, given tlen_bc and a fill sequence.
        """
        # Extract p7 index
        p7 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            ref=settings["ref_p7"],
            targ_dir=settings["tdir_p7"],
            targ_len=settings["tlen_p7"],
            mismatches=settings["mismatches"],
        )

        # If found, set attrs and continue
        if p7[0]:
            self.p7_seq, self.p7_qual = p7[0], np.mean(p7[1])

            # Extract barcode for dual index
            bc = self.match_extract_dual_targ(
                seq=self.seq,
                qual=self.qual,
                ref1=settings["ref_p7"],
                ref2=settings["ref_bc"],
                targ_len=settings["tlen_bc"],
                mismatches=settings["mismatches"],
                fill_seq=settings["fill_seq"],
            )

            if bc[0]:
                self.bc_seq, self.bc_qual = bc[0], np.mean(bc[1])

    @staticmethod
    def Q33conv(qual: str):
        """
        Convert a string of Q33 ascii characters into a list of PHRED scores.
        """
        return [ord(q) - 33 for q in qual]

    @staticmethod
    def revcomp(seq: str):
        """
        Returns the reverse complement of a DNA sequence.
        """
        complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "X": "X"}
        return "".join(complement[base] for base in seq[::-1].upper())

    @staticmethod
    def match_extract_single_targ(
        seq: str, 
        qual: list[int], 
        ref: str, 
        targ_dir: str, 
        targ_len: int,
        mismatches: int,
    ):
        """
        Extract sequence of given length and direction by single target.
        """
        # locate match
        #refmatch = re.search(ref, seq)
        refmatch = regex.search(f"({ref}){{s<={mismatches}}}", seq)

        if refmatch:
            # get positions
            start, end = refmatch.start(), refmatch.end()

            # return appropriate slice
            if targ_dir == "upstream":
                return [seq[start - targ_len : start], qual[start - targ_len : start]]
            elif targ_dir == "downstream":
                return [seq[end : end + targ_len], qual[end : end + targ_len]]

            # or tell user we need a direction
            else:
                raise ValueError("targ_dir param requires 'upstream' or 'downstream'")

        # if no match, return empty
        return ["", []]

    @staticmethod
    def match_extract_dual_targ(
        seq: str,
        qual: list[int],
        ref1: str,
        ref2: str,
        targ_len: int,
        mismatches: int,
        fill_seq: str | None,
    ):
        """
        Extract sequence between ref1 and ref2.
        Pad front of xseq to targ_len with 3' of fill_seq, if fill_seq provided
        """
        # locate matches
        #ref1match, ref2match = re.search(ref1, seq), re.search(ref2, seq)
        ref1match = regex.search(f"({ref1}){{s<={mismatches}}}", seq)
        ref2match = regex.search(f"({ref2}){{s<={mismatches}}}", seq)

        if ref1match and ref2match:
            start, end = ref1match.end(), ref2match.start()

            # if start is before or at end (this catches empty short bcs!)
            if start <= end:
                # get slices
                xseq, xqual = seq[start:end], qual[start:end]

                # pad 5' bc_seq to desired length with 3' of fill_seq
                # (if fill_seq is provided in stock settings, and bc is short)
                if fill_seq and len(xseq) < targ_len:
                    diff = targ_len - len(xseq)
                    xseq = fill_seq[-diff:] + xseq
                    xqual = [30] * diff + xqual

                return [xseq, xqual]

        # if either ref is missed, return empty
        return ["", []]
