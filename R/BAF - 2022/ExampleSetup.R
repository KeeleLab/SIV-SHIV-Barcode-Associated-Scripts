#### Inputs that we need to pass to the function (change frequently) ####

# Path to folder containing fastq file and runinfo excel file
folder <- "~/Documents/BarcodeAnalysisTool/ExampleSetup/"

# Path to the reference FASTA containing known barcode list
# There is a different reference file for each barcoded virus
# reffilepath <- NA # For uncharacterized/new stocks
reffilepath <- "~/Documents/BarcodeAnalysisTool/239M reference 071619.fasta"
#reffilepath <- "~/Documents/BarcodeAnalysisTool/239M2 reference.fasta"
# reffilepath <- "~/Documents/BarcodeAnalysisTool/SIVmac239X Reference Sequences.fasta"

#### Inputs that rarely change unless doing a complicated analysis ####

# Adapter sequence (ours always has NNNN)
adapter = "NNNN"

# Specify whether reads are forward or reverse complemented
fwdrev <- "Reversed" # For 239M and 239M2
#fwdrev <- "Forward" # For X-virus/INT

# Reference sequence reads will be aligned to
ref = "ATGGAAGAAAGACCTCCAGAAAATGAAG" # For 239M and 239M2
#ref = "TAAAAATTTTCGGGTCTATTAC" # For X-virus/INT

# Maximum number of mismatches to reference that are allowed when aligning
mismatches <- 2

# offset is used to specify distance between ref and target. Typically zero
# except for custom analyses
offset <- 0

# Specify whether to extract region upstream or downstream of reference
direction <- "Upstream" # For 239M and 239M2
#direction <- "Downstream" # For X-virus/INT

# How many bases to cut
numbases <- 34 # For 239M and 239M2
# numbases <- 45 # For X-virus/INT

# Whether to use 1/input cutoff
cutoff <- "Yes"

# Whether to include uniques in full analysis or not
full <- FALSE # For normal runs

# Hamming distance to flag
mindist <- 1

# How many reads to analyze at a time (more is faster but uses more memory)
parse <- 5*10^5

# Path to excel file containing Miseq primer sequences
primerseqs <- "~/Documents/BarcodeAnalysisTool/106 primer list.xlsx"

#### Common analysis settings ####

### Normal VPX run settings
# Reference file: 239M reference
# Reverse complemented
# Reference sequence: ATGGAAGAAAGACCTCCAGAAAATGAAG
# Extract upstream of reference
# Extract 34 nucleotides


### Start analysis ####

barcode_analysis(folder, reffilepath, adapter, fwdrev, ref, mismatches, direction, numbases, cutoff, full, mindist, parse, primerseqs)
