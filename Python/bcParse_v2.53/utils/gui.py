# gui.py

"""
Name:       gui.py
Author:     CAG
Version:    1.3
Date:       2025/4/14
"""

# %% Imports
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# %% Classes


class ArgSelectionGUI:

    def __init__(self, stock_list: list[str]):
        # Create the main window
        self.root = tk.Tk()
        self.root.title("File and Stock Selection")

        # Create a main frame - autoexpands vs root.geometry
        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        # Store stock list for dropdown
        self.stock_list = stock_list

        # Dictionary to store selections
        self.args = {}

        # Dictionary to store widgets labels
        self.labels = {}

        self.out_label = tk.Label(self.root)
        self.stock_var = tk.StringVar(self.root)

        # Define file types and their corresponding dialog options
        self.file_types = {
            "sample.fq": ("Select sample fastq file", "Fastq files", ("*.*")),
            "runinfo.xlsx": ("Select Runinfo Excel file", "Excel files", "*.*"),
        }

        # Initialize BooleanVars with defaults
        self.filt_mat_ac_var = tk.BooleanVar(value=True)
        self.mask_var = tk.BooleanVar(value=True)
        self.legacy_format_var = tk.BooleanVar(value=True)

        # Initialize IntVars with defaults
        self.mean_qual_var = tk.IntVar(value=30)
        self.mismatches_var = tk.IntVar(value=1)
        self.mask_qual_var = tk.IntVar(value=20)

        self.setup_ui()

    def browse_file(self, file_type):
        """
        Define action for file selection buttons
        """
        # Open file dialog and update file path
        title, file_desc, file_ext = self.file_types[file_type]

        path = filedialog.askopenfilename(
            title=title, filetypes=[(file_desc, file_ext)], initialdir=os.getcwd()
        )

        if path:
            self.args[file_type] = path
            self.labels[file_type].config(text=f"{file_type}: {path}")

    def browse_output(self):
        """
        Define action for output directory button
        """
        # Open directory dialog for output selection
        path = filedialog.askdirectory(
            title="Select Output Directory", initialdir=os.getcwd()
        )

        if path:
            self.args["out_path"] = path
            self.out_label.config(text=f"output: {path}")

    def on_submit(self):
        """
        Define action for submit button
        """
        # Check if all required files have been selected
        required_fields = ["sample.fq", "runinfo.xlsx", "out_path"]

        # Collect missing fields
        missing_fields = [field for field in required_fields if field not in self.args]

        # Return message if missing.
        if missing_fields:
            missing_str = ", ".join(missing_fields)
            messagebox.showerror(
                "Missing Selections", f"Please select the following: {missing_str}"
            )
            return

        # Add stock selection & other params to args if all are present
        self.args["stock"] = self.stock_var.get()

        self.args["mask"] = self.mask_var.get()
        self.args["filt_mat_ac"] = self.filt_mat_ac_var.get()
        self.args["legacy_format"] = self.legacy_format_var.get()

        self.args["mean_qual"] = self.mean_qual_var.get()
        self.args["mismatches"] = self.mismatches_var.get()
        self.args["mask_qual"] = self.mask_qual_var.get()

        # All required selections are made, close the window
        self.root.quit()

    def setup_ui(self):
        """
        Build the GUI
        """
        # Create a subframe for required selections
        self.req_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.req_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)  # Adjust packing as needed

        # Add a label to the subframe
        req_label = tk.Label(self.req_frame, text="Required selections")
        req_label.pack(padx=2, pady=2, anchor="w")

        # Create selection buttons -
        # For files,
        for _, (file_type, _) in enumerate(self.file_types.items()):
            frame = tk.Frame(self.req_frame)
            frame.pack(fill=tk.X, expand=True, pady=2)

            tk.Button(
                frame,
                text=f"Browse {file_type}",
                command=lambda ft=file_type: self.browse_file(ft),
            ).pack(side=tk.LEFT, padx=5)

            self.labels[file_type] = tk.Label(frame, text=f"{file_type}: Not selected")

            self.labels[file_type].pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # For output directory,
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Button(frame, text="Browse output", command=self.browse_output).pack(
            side=tk.LEFT, padx=5
        )

        self.out_label = tk.Label(frame, text="output: Not selected")

        self.out_label.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # For stock selection drop-down,
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Select stock:").pack(side=tk.LEFT, padx=5)

        # Variable to hold the stock selection, default 239M
        self.stock_var.set("239M_dual-index")

        stock_menu = ttk.Combobox(
            frame, textvariable=self.stock_var, values=self.stock_list
        )

        stock_menu.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=False)

        # Create a subframe for optional selections
        self.opt_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.opt_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)  # Adjust packing as needed

        # Add a label to the subframe
        opt_label = tk.Label(self.opt_frame, text="Optional - Fastq QC")
        opt_label.pack(padx=2, pady=2, anchor="w")

        # Add checkbox for mask
        frame = tk.Frame(self.opt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.mask_checkbox = tk.Checkbutton(
            frame, text="Mask low-q bases to collapse/drop Ns", variable=self.mask_var
        )
        self.mask_checkbox.pack(side=tk.LEFT, padx=2)

        # Add spinbox for mask_qual
        frame = tk.Frame(self.opt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Mask PHRED below:").pack(side=tk.LEFT, padx=2)

        self.mask_qual_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=40,
            textvariable=self.mask_qual_var,
            width=5
        )
        self.mask_qual_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

        # Add spinbox for mean_qual
        frame = tk.Frame(self.opt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Drop PHRED below:").pack(side=tk.LEFT, padx=2)

        self.mean_qual_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=40,
            textvariable=self.mean_qual_var,
            width=5
        )
        self.mean_qual_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

         # Add spinbox for mismatches
        frame = tk.Frame(self.opt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Allow mismatches?").pack(side=tk.LEFT, padx=2)

        self.mismatches_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=5,
            textvariable=self.mismatches_var,
            width=5
        )
        self.mismatches_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

        # Create a subframe for optional selections
        self.fmt_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.fmt_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)  # Adjust packing as needed

        # Add a label to the subframe
        fmt_label = tk.Label(self.fmt_frame, text="Optional - Excel format")
        fmt_label.pack(padx=2, pady=2, anchor="w")

        # Add checkbox for legacy_format
        frame = tk.Frame(self.fmt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.legacy_format_checkbox = tk.Checkbutton(
            frame,
            text="Use legacy format",
            variable=self.legacy_format_var,
        )
        self.legacy_format_checkbox.pack(side=tk.LEFT, padx=2)

        # Add checkbox for filt_mat_ac
        frame = tk.Frame(self.fmt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.filt_mat_ac_checkbox = tk.Checkbutton(
            frame,
            text="Filter matrix above-cutoff",
            variable=self.filt_mat_ac_var,
        )
        self.filt_mat_ac_checkbox.pack(side=tk.LEFT, padx=2)

        # And for Submit.
        submit_button = tk.Button(
            self.main_frame, text="Submit", command=self.on_submit
        )
        submit_button.pack(pady=10, anchor="w")

    def run(self) -> tuple[str, str, str, str, bool, bool, bool]:
        """
        Run the GUI and return a list of args
        """
        self.root.mainloop()  # Wait for user interaction

        try:
            self.root.destroy()  # Attempt to safely destroy the root window
        except tk.TclError:
            pass  # Window was already destroyed, ignore the error

        sample_path = str(self.args.get("sample.fq"))
        runinfo_path = str(self.args.get("runinfo.xlsx"))
        out_path = str(self.args.get("out_path"))
        stock = str(self.args.get("stock"))
        mean_qual = int(self.args.get("mean_qual"))
        mismatches = int(self.args.get("mismatches"))
        mask = bool(self.args.get("mask"))
        mask_qual = int(self.args.get("mask_qual"))
        filt_mat_ac = bool(self.args.get("filt_mat_ac"))
        legacy_format = bool(self.args.get("legacy_format"))

        return (
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
        )
