#### Needed Packages ####################################################################################################

# Before running anything, load each of these packages

# Check for needed packages and install if necessary, may take a few minutes
if (!"shiny" %in% rownames(installed.packages())){
  install.packages('shiny',dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"rhandsontable" %in% rownames(installed.packages())){
  install.packages('rhandsontable', dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"stringr" %in% rownames(installed.packages())){
  install.packages('stringr', dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"ShortRead" %in% rownames(installed.packages())){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ShortRead")}
if (!"shinyFiles"%in%rownames(installed.packages())){
  install.packages("shinyFiles",dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"shinythemes"%in%rownames(installed.packages())){
  install.packages("shinythemes",dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"stringdist"%in%rownames(installed.packages())){
  install.packages("stringdist",dependencies=TRUE, repos='http://cran.rstudio.com/')}
if(!"readxl"%in%rownames(installed.packages())){
  install.packages("readxl", dependencies = TRUE, repos = 'http://cran.rstudio.com/')}
if(!"openxlsx"%in%rownames(installed.packages())){
  install.packages("openxlsx", dependencies = TRUE, repos = "http://cran.rstudio.com/")}
if(!"anytime" %in% rownames(installed.packages())){
  install.packages("anytime", dependencies = TRUE, repos = "http://cran.rstudio.com/")
}

library(shinythemes)
library(shiny)
library(ShortRead)
library(stringr)
library(rhandsontable)
library(shinyFiles)
library(stringdist)
library(readxl)
library(openxlsx)
library(anytime)

  #### 05/05/2022 - Barcode App w/ Short Seqs Flagged, Compatible with all current analysis runs. - ANC
  
  # Changes across versions:
  #### 05/05/2022 -
  # 1. Commented out lines of code that did not add reads to the analysis file if seqs.Count
  # only had a single row. Affects primarily low input data with a low number of barcodes reported.
  # 2. Corrected an issue with capitalization of primer indices. VPX/vpx/vPx and INT,int,inT should
  #     now all be valid inputs on the run info sheet. 
  # 3. Changed id() to ShortRead::id() to avoid package conflicts
  # 02/08/2022
  # 1. Changes to primer index gsub to incorporate errors made when creating new primers for Nef samples.
  #     Needed to sub TAGAA without subbing from CTAGAA
  # 11/30/2021
  # 1. Fixed an issue where the folder name would not properly be added to the name of the output file.
  
barcode_analysis <- function(folder, reffilepath, adapter = "NNNN", fwdrev = "Reversed", ref = "ATGGAAGAAAGACCTCCAGAAAATGAAG", mismatches = 2, direction = "Upstream", numbases = 34, cutoff = "Yes", full = F, mindist = 1, parse = 5*10^5, primerseqs){
  
  if(!dir.exists(folder)){
    stop(paste0("Folder ", folder, " does not exist"))
  }else if(!is.na(reffilepath)){
    if(!file.exists(reffilepath)){
      stop(paste0("Reference file ", reffilepath, " does not exist"))
    }
  }else if(!file.exists(primerseqs)){
    stop(paste0("Miseq primer sequences in file ", primerseqs, " cannot be found"))
  }
  
  ## Grab parent directory from reffilepath to find other barcode stocks
  if(!is.na(reffilepath)){
    reffilepath_directory <- unlist(strsplit(reffilepath, split="/"))
    reffilepath_directory <- paste0(reffilepath_directory[1:length(reffilepath_directory)-1], collapse = "/")
    barcodeStock <- unlist(strsplit(reffilepath, split="/"))[length(unlist(strsplit(reffilepath, split="/")))]
    VPXVPRstocks <- c("239M reference 071619.fasta", "239M2 reference.fasta", "NIR  barcode list fasta.fasta", "OptM Reference.fasta", "PolH41Y I50V and 239M combined barcode list.fasta")
    if(barcodeStock %in% VPXVPRstocks){
      ContamCheck <- TRUE
      correctIndex <- which(barcodeStock == VPXVPRstocks)
      barcodeStock <- c("239M", "293M2", "NIRM", "optM", "Pol")
      incorrectStocks <- barcodeStock[-correctIndex]
      barcodeStock <- barcodeStock[correctIndex]
      if(barcodeStock == "239M"){
        reffile <- readFasta(reffilepath_directory, "239M reference 071619.fasta")
      }else if(barcodeStock == "293M2"){
        reffile <- readFasta(reffilepath_directory, "239M2 reference.fasta")
      }else if(barcodeStock == "NIRM"){
        reffile <- readFasta(reffilepath_directory, "NIR  barcode list fasta.fasta")
      }else if(barcodeStock == "optM"){
        reffile <- readFasta(reffilepath_directory, "OptM Reference.fasta")
      }else if(barcodeStock == "Pol"){
        reffile <- readFasta(reffilepath_directory, "PolH41Y I50V and 239M combined barcode list.fasta")
      }
      combinedVPXStocks <- readFasta(paste0(reffilepath_directory, "/CombinedVPXStocks.fasta"))
      combinedSeqs <- as.character(sread(combinedVPXStocks))
      combinedNames <- as.character(ShortRead::id(combinedVPXStocks))
      reffileseq <- as.character(sread(reffile))
      reffilename <- as.character(ShortRead::id(reffile))
      combinedNames <- combinedNames[-which(combinedSeqs %in% reffileseq)]
      combinedSeqs <- combinedSeqs[-which(combinedSeqs %in% reffileseq)]
    }else{
      print("Stock not found in VPX list, will continue as non-standard stock analysis")
      ContamCheck <- FALSE
      reffile <- readFasta(reffilepath)
      reffileseq <- as.character(sread(reffile))
      reffilename <- as.character(ShortRead::id(reffile))
    }
  }else{
    ContamCheck <- FALSE
    reffileseq <- NA
    reffilename <- NA
  }
  version <- "04/18/2022 - Across Stock Validation version"
  
  # If the user specified less that 10,000 reads at a time, restore to 500,000 default
  if (is.na(as.numeric(parse))|(as.numeric(parse) < 10^4)){parse = 5*10^5}
  
  # Table on Tab 2 listing all primers
  
  # First row: Primer index names, 0 is for if no primer was used, 1:40 are our usual sequencing primers
  # The rest are for filling in if a primer was used that isn't one of these options
  primer <- paste("P5.", 0:399, sep = "")
  
  # Second row: Primer sequences
  # Read in sequences from file
  primerseqs <- suppressMessages(read_excel(primerseqs, col_names = FALSE))
  primerseqs <- as.data.frame(primerseqs)
  primerseqs[,2] <- gsub("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNN|CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCNNNN","",
                         unlist(primerseqs[,2]))
  primerseqs[,2] <- gsub("CCAGAACCTCCACTACCCATTCATC|AGAAGCGGAGACAGCGACGAAG|CGCCCTTACTGCCTTCACTCAGC|CTCTGATAATCTGCATAGCCGCTTG|CAGCGAGTTTCCTTCTTGTCAGCC|CCTCTACTGTTCTGTTTAGCCATTCG|CTTGCTTTACAGCGGGAGAAGTG|GTGACAGAGGATTTGCTGCACC|GAACTCATTAGAATCCTCCAACGAGCG|GTGACTGTATTATTCCAGGCATCAAAGC|CCCTCTTATTTCCAGCAGACCCATA|GATTGCGAGTATCCATCTTCCACCTCTC|CAGCGAGTTTCCTTCTTGTCAGCC|GAGTTTGGAAGCAAGTCAGGCCTG|GGATGGGCATGGTGGACCTG|GCGACCCTACAGAGGATTCGAGAAG|ACGCGACCCTACAGAGGATTCGAGAAG","", unlist(primerseqs[,2]))
  primerseqs[,2] <- toupper(unlist(primerseqs[,2]))
  
  ##ignore this for distributed version, these lines were written for a custom analysis
  
  # nefPrimers <- grep("TAGAACCTCTCCCCAAGGGTC", primerseqs[,2])
  # misalignedIndex <- grep("GTACTGAC|TGAATACC|GCCATGAC", primerseqs[,2])
  # bothNefAndMisaligned <- nefPrimers[nefPrimers %in% misalignedIndex]
  # primerseqs[,2][bothNefAndMisaligned] <- gsub("TAGAACCTCTCCCCAAGGGTC", "", primerseqs[,2][bothNefAndMisaligned])
  # primerseqs[nefPrimers,2] <- gsub("CTAGAACCTCTCCCCAAGGGTC|TAGAACCTCTCCCCAAGGGTC", "", primerseqs[nefPrimers,2])
  
  primerseqs <- primerseqs[order(as.numeric(gsub("VPX.P5.|VPU.P5.|PBS.P5.|Gag.P5.|Nef.P5.|VifRL8.P5.|VifRL9.P7.|PR.P5.|CTL.amp1.P5.|CTL.amp1.R2.P5.|CTL.amp2.P5.|CTL.amp2.R2.P5.|Nef.RL9.Rev.P5.|Nef.RL9.R2.P5.|INT.P5.|HW8.P5.|LY10.P5.|Nef.escapes.P5.", "", unlist(primerseqs[,1])))),]
  # 40 multiplexed sequencing adapters (primer indexes), first one is empty for already-demultiplexed data, rest can be filled in with custom sequences
  fillto = 399-nrow(primerseqs)
  primerseqs <- c("", unlist(primerseqs[,2]), rep("", fillto))
  
  # Third row: TRUE/FALSE, table will have checkboxes for user to select primers
  Include = rep(FALSE, 400)
  
  #### Match up  runinfo and fastq files ####
  # For batch version, read all files in a folder as input instead of single fastq file
  files <- dir(path = folder, recursive = TRUE)
  # If any files are open the temp file will also show up...
  files <- files[!grepl("[~$]|Analysis|depreciated", files, ignore.case=T)]
  
  # List all the files in this folder
  #print(files)
  # The xlsx files should be our runinfo files
  batchinfo <- data.frame(runinfo = files[grepl(".xlsx", files)], `Run Name` = NA, fastq = NA,
                          stringsAsFactors = FALSE)
  
  # How do we match up the runinfo files with their fastq files?
  # We should be able to find the line in the runinfo file that contains the fastq filename
  batchinfo$Run.Name <- unlist(lapply(batchinfo$runinfo, function(x){
    info <- suppressMessages(read_excel(paste0(folder, "/", x)))
    column <- grep("Run Name", info)
    row <- grep("Run Name", as.character(unlist(info[,column])))
    as.character(unlist(info[,column]))[row + 1]
  }))
  
  # Can't have colons in excel sheet name so check for that
  if(TRUE %in% grepl(":", batchinfo$runinfo)){
    stop("Can't have ':' or '/' in runinfo file name")
  }
  
  # And then we can match these up to the fastq files
  for (i in files[grepl(".fastq|fq", files)]){
    match <- grepl(gsub("[_].+", "", basename(i)), batchinfo$Run.Name, ignore.case = TRUE)
    if(sum(match) == 0){
      # If we can't identify the file print an error
      stop("unable to match fastq file ", basename(i), " to runinfo")
    }else if(sum(match > 1)){
      stop("multiple fastq files named ", basename(i))
    }
    batchinfo$fastq[match] <- i
  }
  
  # Print out what files were matched up just to check
  for(i in 1:nrow(batchinfo)){
    cat(paste0('Analyzing run "', batchinfo[i,2], '" with runinfo file "', batchinfo[i,1],
                 '" and fastq file "', gsub(".+[/]", "", batchinfo[i, 3])), '"\n')
  }
  
  # Check each run for proper formatting before starting the full analysis
  for(file in 1:nrow(batchinfo)){
    
    # Load runinfo file
    if(nrow(batchinfo) == 1){
      runinfo <- suppressMessages(read_excel(paste0(folder, "/", as.character(batchinfo$runinfo))))
    }else{
      runinfo <- suppressMessages(read_excel(paste0(folder, "/", as.character(unlist(batchinfo$runinfo))[file])))
    }
    
    # And figure out which primers we have
    if("Barcodes" %in% colnames(runinfo)){
      Include <- unlist(lapply(primer, function(x){sum(grepl(paste0(x, "$"), runinfo$Barcodes, ignore.case = TRUE)) > 0}))
    }else{
      stop(paste0("Unable to find primers for file ", batchinfo$runinfo[file]))
    }
    
    # As well as what their sequencing inputs were
    Input = rep(0, 400)
    Input[Include] <- unlist(lapply(primer[Include], function(x){
      runinfo$`Input TOTAL PER BARCODE`[which(grepl(paste0(x, "$"), runinfo$Barcodes))]
    }))
    
    # Dataframe holding this information for the table
    # Primer number, sequence, whether it was used, and input value
    primer.list <- data.frame(primer, primerseqs, Include, Input, stringsAsFactors = FALSE)
    
    # Column names, will also be headers in the table
    colnames(primer.list) = c("5' Index Name", "Index Sequence (blank if none)", "Include?",
                              "Sequencer Template Input (0 if unknown or unused)")
    
    # Row names to display in the table
    rownames(primer.list) <- 0:399
    
    # Other potential checks:
    # Is sequencing info in correct format (numerical not character)
    # Are runinfo column headers correct (certain headers need to match exactly to what code expects)
    
  }
  
  # This is where we will save our excel tabs/sheets
  list_of_datasets <- list()
  # And this will be the full matrix version for Taina
  list_of_datasets_full <- list()
  
  #### Loop through each fastq file ####
  for(file in 1:nrow(batchinfo)){
    
    #### Set up miseq run ####
    
    # Load runinfo file
    if(nrow(batchinfo) == 1){
      runinfo <- suppressMessages(read_excel(paste0(folder, "/", as.character(batchinfo$runinfo))))
    }else{
      runinfo <- suppressMessages(read_excel(paste0(folder, "/", as.character(unlist(batchinfo$runinfo))[file])))
    }
    
    # And figure out which primers we have
    if("Barcodes" %in% colnames(runinfo)){
      Include <- unlist(lapply(primer, function(x){sum(grepl(paste0(x, "$"), runinfo$Barcodes, ignore.case = TRUE)) > 0}))
    }else{
      stop(paste0("Unable to find primers for file ", batchinfo$runinfo[file]))
    }
    
    # As well as what their sequencing inputs were
    Input = rep(0, 400)
    Input[Include] <- unlist(lapply(primer[Include], function(x){
      runinfo$`Input TOTAL PER BARCODE`[which(grepl(paste0(x, "$"), runinfo$Barcodes))]
    }))
    
    # Dataframe holding this information for the table
    # Primer number, sequence, whether it was used, and input value
    primer.list <- data.frame(primer, primerseqs, Include, Input, stringsAsFactors = FALSE)
    
    # Column names, will also be headers in the table
    colnames(primer.list) = c("5' Index Name", "Index Sequence (blank if none)", "Include?",
                              "Sequencer Template Input (0 if unknown or unused)")
    
    # Row names to display in the table
    rownames(primer.list) <- 0:399
    if (fwdrev == "Forward") {
      primer_prefix = "INT."
    } else if (fwdrev == "Reversed") {
      primer_prefix = "VPX."
    } else {
      stop("fwdrev not set")
    }
    runinfo$Barcodes[!is.na(runinfo$Animal)] <- toupper(runinfo$Barcodes[!is.na(runinfo$Animal)])
    P5barcodes <- regexec("[0-9]+$", runinfo$Barcodes)
    P5barcodes <- as.numeric(regmatches(runinfo$Barcodes, P5barcodes))
    if(sum(P5barcodes[!is.na(P5barcodes)] != P5barcodes[order(P5barcodes)][!is.na(P5barcodes)]) == length(P5barcodes[!is.na(P5barcodes)])){
      probPrimers <- runinfo$Barcodes[which(P5barcodes != P5barcodes[order(P5barcodes)])]
      stop(paste("Implement formatting for primers listed out of order, or atypical primer prefix used in runinfo file (VPX or Int?) See following:", paste(probPrimers, collapse = " ")))
    }
    
    #### Analysis Setup####
    
    #for(qscore in 15:30){
    # 2/17/20 Q30 trimming added
    qscore <- 30
      
    # Settings and inputs from user
    
    # Table with adapters (primers) that the user selected on Tab 2 (or default table)
    #primer.list.cleaned = DF()
    primer.list.cleaned <- primer.list
    
    # Path to the file we are analyzing
    #path.to.file = as.character(parseFilePaths(volumes, input$runfile)$datapath)
    if(nrow(batchinfo) == 1){
      path.to.file <- paste0(folder, "/", batchinfo$fastq)
    }else{
      path.to.file <- paste0(folder, "/", batchinfo$fastq[file])
    }
    
    # Percent contamination (index hopping) to flag (Tab 3)
    # Can determine whether to have the index hopping flag
    # By seeing if there are multiple animals or not
    if(length(unique(runinfo$Animal[!is.na(runinfo$Animal)])) > 1 & !is.na(reffilepath)){
      contam = 0.2
    }else{
      contam = 100
    }
    
    # Set up data stream to pull reads from it in chunks of the specified size
    stream = FastqStreamer(path.to.file, parse)
    # Close connection when we're done so it doesn't throw an error
    on.exit(close(stream))  
  
    # Set up tables to hold our data for each primer
    
    # If the user selected primers, only use those specified primers
    if (sum(primer.list.cleaned[,3]) > 0){
      primer.list.cleaned = primer.list.cleaned[primer.list.cleaned[,3],]
    }else{
      
      # If the user didn't specify any primers
      # Then we'll run all the data as one sample so we only need one row in the table
      primer.list.cleaned = primer.list.cleaned[1,];
      primer.list.cleaned[1,2] = ""
      
    }
    
    # Convert barcode (indexing) sequences into a character vector
    primer.list.cleaned[,2] = as.character(primer.list.cleaned[,2])
    # Note if no sequencing input is given for a primer
    primer.list.cleaned[is.na(primer.list.cleaned[,4]), 4] = Inf
    primer.list.cleaned[primer.list.cleaned[,4] == 0,4] = Inf
    
    # Create a big holder matrix that we will put all of our analysis results into
    barcodeunique.hold <- as.data.frame(matrix(NA, nrow = 1*10^4, ncol = (nrow(primer.list.cleaned)*8)))
    # Second row of this matrix will have the sequencing input, the number of sequences found per primer/in this sample
    # The number of barcodes found in this sample, and what proportion of reads were unique/previously unknown barcodes
    # And what the index hopping threshold was if used
    barcodeunique.hold[2,] = rep(c("Sequencer Input", "Total # Sequences", "# of Barcodes (Excluding Uniques)",
                                   "Proportion Uniques", "Proportion Contamination", NA, NA, NA),
                                 ncol(barcodeunique.hold)/8) 
    # Fifth row of the matrix lists the barcode name, sequence, counts, proportion out of the sample
    # As well as whether it was flagged for hamming distance and if so to what other barcode
    # And whether it was flagged for possible index hopping
    barcodeunique.hold[5,] = rep(c("Barcode", "Sequence", "Counts", "Proportion", "Hamming Dist to Another Sequence",
                                   "Index hopping?", "Short barcode?", NA), ncol(barcodeunique.hold)/8)
    
    # Create a matrix to contain summary statistics about each primer/sample
    sequencestats = as.data.frame(matrix(0, nrow = nrow(primer.list.cleaned), ncol=6))
    # This matrix will have info about how many sequences were found with each primer,
    # how many aligned to the reference sequence and had a barcode extracted,
    # how many of those extracted sequences matched to a known barcode from our reference file
    # the total number of known barcodes found for each primer, and what the sequencing input was
    colnames(sequencestats) = c("Primer", "# of 5' Indexing Sequences", "Total Sequences Extracted per Primer",
                                "# of Sequences Matching Known Barcode", "# of Known Barcodes", "Sequencing Input")
    # Use the primers specified by the user, or the default
    sequencestats[,1] = primer.list.cleaned[,1]
    
    # #### Analysis Run ###############################################################################################
    
    # Count to keep track of how many sequences have been analyzed
    sum = 0
  
    # Pull first batch of reads
    while(length(fq <- yield(stream))){
      
      # Pull reads from fastq class object
      file.seq = as.character(sread(fq))
      
      # 2/17/20
      # How many reads pass strict Q>30 filter?
      # Going to include quality scores at least until we get to the barcode extraction stage
      
      # Pull Q scores as well
      qual.seq <- as.character(quality(fq)@quality)
      
      # Remove sequencing adapters
      
      # If the adapter is only Ns, take out the length of the adapter from each read
      if (sum(str_detect(adapter, LETTERS[-14])) == 0){
        file.seq <- str_sub(file.seq, nchar(adapter) + 1, nchar(file.seq))
        # And remove these bases from qual score as well
        qual.seq <- str_sub(qual.seq, nchar(adapter) + 1, nchar(qual.seq))
        
      }else{
        # If the adapter has letters other than N
        # Cut out where we expect the adapter to be
        check.adapter <- str_sub(file.seq, 1, nchar(adapter))
        # Keep the reads only if this piece matches the adapter we want
        file.seq <- file.seq[check.adapter == adapter]
        qual.seq <- qual.seq[check.adapter == adapter]
        # Remove the adapter from the reads we're keeping
        file.seq <- str_sub(file.seq, nchar(adapter) + 1, nchar(file.seq))
        qual.seq <- str_sub(qual.seq, nchar(adapter) + 1, nchar(qual.seq))
      }
      
      # Reallocate memory
      # Fastq files contain the read names, sequences, and quality
      # We're freeing up memory used for the metadata that we aren't using in our analysis
      gc() 
      
      ## Demultiplexing ####
      # Cut out the section where we expect the primer sequence to be
      if(FALSE %in% (nchar(primer.list.cleaned[,2]) == nchar(primer.list.cleaned[1,2]))){
        stop("Primer sequences have variable length")
      }else{
        file.seq.4bpTrim <- str_sub(file.seq, 0, nchar(primer.list.cleaned[1,2]))
      }
      
      # Analyze reads one primer/sample at a time
      for (i in 1:nrow(primer.list.cleaned)){
        # Find reads where the primer section matches the sample we're considering
        index = file.seq.4bpTrim == primer.list.cleaned[i,2]
        
        # If less than 0.01% of these reads match the primer, move on to the next one
        # Unless we have the sequencing input and know that there should be reads for this sample
        # In which case continue so long as there are any reads in there at all
        if (((sum(index) < (0.00001*length(file.seq.4bpTrim)/nrow(primer.list.cleaned)))&(primer.list.cleaned[i,4] == Inf)) | (sum(index) == 0 & primer.list.cleaned[i,4] > 0)){
          
          # But still count how many sequences were found in this sample
          sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
          
          next
        }
        
        # Cut the reads containing the primer from right after the primer sequence to the end of the read
        file.seq.demux <- str_sub(file.seq[index], nchar(primer.list.cleaned[i,2]) + 1, nchar(file.seq[index]))
        qual.seq.demux <- str_sub(qual.seq[index], nchar(primer.list.cleaned[i,2]) + 1, nchar(qual.seq[index]))
        
        # Once this is done, remove already demultiplexed reads from file.seq list
        # Matching primers takes time and we already know what these are
        file.seq <- file.seq[!index]
        qual.seq <- qual.seq[!index]
        file.seq.4bpTrim <- file.seq.4bpTrim[!index]
        
        # If these are the reverse strand, reverse complement the reads before we look for the reference sequence
        if (fwdrev == "Reversed"){
          file.seq.demux <- as.character(reverseComplement(DNAStringSet(file.seq.demux)))
          # And use UTF-8 magic to reverse the order of the Q-scores too (this is faster than split and paste)
          qual.seq.demux <- lapply(qual.seq.demux, function(x){intToUtf8(rev(utf8ToInt(x)))})
          
        }
        
        ## Barcode extraction ####
        ## For when the barcode is upstream of the offset
        if (direction == "Upstream") {
          # Match the reads to the reference sequence, allowing the selected number of mismatches
          index.downstream.num = aregexec(ref, file.seq.demux, max.distance = mismatches)
          # Unlist the matches to just get the position of the reference in the read
          index.downstream.num = unlist(index.downstream.num)
          noMatches <- which(index.downstream.num == -1)                        # find where aregexec is equal to -1, indicating no match of reference in seq
          sapply(noMatches, function(x){
            index.downstream.num[[x]] <<- -10000                                # set index.downstream.num to large negative number for these so it gets filtered out in next step
          })
          # We're only interested in sequences where the reference is found far enough from the beginning of the read
          # For there to be space for the whole barcode/region of interest to be there
          file.seq.demux.noNA <- file.seq.demux[index.downstream.num + offset >= numbases + 1]
          qual.seq.demux.noNA <- qual.seq.demux[index.downstream.num + offset >= numbases + 1]
          # Remove the reads we don't want from the index list as well
          index.downstream.num.noNA = index.downstream.num[index.downstream.num + offset >= numbases + 1]
          
          # If less than 2 reads meet these criteria, move on to the next primer
          if (length(index.downstream.num.noNA) < 2) {
            
            # But still count how many sequences were found in this sample
            sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
            
            next
          }
          if (offset != 0){
          # Cut out the section where we expect our barcode to be (default is the 34 bases before the reference sequence)
          seqs = str_sub(
            file.seq.demux.noNA,
            index.downstream.num.noNA - numbases + offset,
            index.downstream.num.noNA - 1 + offset
          )
          quals <- str_sub(
            qual.seq.demux.noNA,
            index.downstream.num.noNA - numbases + offset,
            index.downstream.num.noNA - 1 + offset
          )
          } else if (offset == 0){
              seqs = str_sub(
                file.seq.demux.noNA,
                index.downstream.num.noNA - numbases,
                index.downstream.num.noNA - 1
              )
              quals <- str_sub(
                qual.seq.demux.noNA,
                index.downstream.num.noNA - numbases,
                index.downstream.num.noNA - 1
              )
            }
          # Convert these sequences to a dataframe
          seqs = as.data.frame(DNAStringSet(seqs))
          
          # And calculate the average per base read quality within the barcode region
          means <- unlist(lapply(quals, function(x){
            length(which((utf8ToInt(x) - 33) < 30)) }))
          
          # Remove barcodes that had scores less than Q30
          seqs <- seqs[(means <= 1),]
          seqs <- as.data.frame(seqs)
          
        }
        
        # For when the barcode is downstream of the reference
        if (direction == "Downstream") {
          
          # And the reads
          #file.seq.demux = as.character(reverseComplement(DNAStringSet(file.seq.demux)))
          # And their quality scores
          #qual.seq.demux <- lapply(qual.seq.demux, function(x){intToUtf8(rev(utf8ToInt(x)))})
          
          # Match to the reference sequence
          index.downstream.num = aregexec(ref, file.seq.demux, max.distance = mismatches)
          # Unlist to get the starting position of each match
          index.downstream.num = unlist(index.downstream.num)
          seqLength <- nchar(file.seq[2])
          noMatches <- which(index.downstream.num == -1)                        # finds where aregex = -1 indicating no match. 
          sapply(noMatches, function(x){
            index.downstream.num[[x]] <<- 10000                                 # set -1 to a large number so it gets filtered out in the following step
          })
          # Only interested in sequences where the reference is far enough from the start for the barcode to be there
          file.seq.demux.noNA = file.seq.demux[seqLength - index.downstream.num >= numbases + 1]
          qual.seq.demux.noNA <- qual.seq.demux[seqLength - index.downstream.num >= numbases + 1]
          # And only want the indices for those reads as well
          index.downstream.num.noNA = index.downstream.num[seqLength - index.downstream.num >= numbases + 1]
          
          # If less than two reads meet these conditions, move on to the next primer
          if (length(index.downstream.num.noNA) < 2) {
            
            print("No reads in this sample matched reference sequence (checkpoint 2)")
            print(i)
            
            # But still count how many sequences were found in this sample
            sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
            
            next
          }
          
          if (offset != 0){                                                     # offset is a position relative to start of reference. Negative signifies upstream of ref, positive is downstream of ref
          # Extract the section where we expect the barcode to be, defaults to the 34 bases before the reference
          seqs = str_sub(
            file.seq.demux.noNA,
            index.downstream.num.noNA + offset,
            index.downstream.num.noNA + offset + numbases - 1
          )
            quals <- str_sub(
            qual.seq.demux.noNA,
            index.downstream.num.noNA + offset,
            index.downstream.num.noNA + offset + numbases - 1
            )
          } else if (offset == 0) {
            seqs = str_sub(
              file.seq.demux.noNA,
              index.downstream.num.noNA + nchar(ref),
              index.downstream.num.noNA + nchar(ref) + numbases - 1
          )
            quals <- str_sub(
              qual.seq.demux.noNA,
              index.downstream.num.noNA + nchar(ref),
              index.downstream.num.noNA + nchar(ref) + numbases - 1
            )
          }
          
          # Convert these sequences to a dataframe
          seqs = as.data.frame(DNAStringSet(seqs))
          # And calculate the average per base read quality within the barcode region
          means <- unlist(lapply(quals, function(x){
            length(which((utf8ToInt(x) - 33) < 30)) }))
          
          # Remove barcodes that had scores less than Q30
          seqs <- seqs[(means <= 1),]
          seqs <- as.data.frame(seqs)
          
        }
        
        # We want to keep track of stats about the data for our sequence stats table
        
        # Make a table of barcodes and how many reads they were in, convert to data frame
        seqs.Count = as.data.frame(table(seqs))
        # Memory allocation, don't need all of the sequences now that they're counted
        rm(seqs); gc()
        
        # # If less than two barcodes found in this primer, move on to the next primer
        # if (nrow(seqs.Count) < 2){
        # 
        #   # But still count how many sequences were found in this sample
        #   sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
        # 
        #   print("Didn't find any reads in this sample (checkpoint 3)")
        #   print(i)
        # 
        #   # next
        # }
        # Convert from factor, sequences to character and counts to numeric
        if (length(seqs.Count) < 2){
          next
        }
        seqs.Count[,1] = as.character(seqs.Count[,1])
        seqs.Count[,2] = as.numeric(as.character(seqs.Count[,2]))
        
        # See what we've already found for this primer in previous batches of reads
        alreadyleft = barcodeunique.hold[6:(nrow(barcodeunique.hold)), (8*i-6):(8*i-5)]
        # Remove any rows containing NAs from the data
        alreadyleft = alreadyleft[!is.na(alreadyleft[,1]),]
        # Convert counts to numeric
        alreadyleft[,2] = as.numeric(alreadyleft[,2])
        # Add new counts from the current chunk of reads to whatever has already been found
        newSeqs.Count = data.frame(tapply(c(seqs.Count[,2], alreadyleft[,2]), c(seqs.Count[,1], alreadyleft[,1]),
                                          FUN = sum))
        # Make sequences into their own column instead of just row names
        newSeqs.Count[,2] = row.names(newSeqs.Count)
        # Reorder columns
        newSeqs.Count = newSeqs.Count[,c(2,1)]
        # And update the table of sequences and counts
        seqs.Count = newSeqs.Count
        # Memory allocation, don't need the temp one anymore
        rm(newSeqs.Count)
        gc()
        
        # If we haven't found any barcodes yet in this primer, move on to the next one
        # if (nrow(seqs.Count) < 2){
  
          # # But still count how many sequences were found in this sample
          # sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
          # 
          # print("Didn't find any reads in this sample (checkpoint 4)")
          # print(i)
  
          # next
        # }
        
        # Otherwise, organize the data for this primer
        
        # Match counted barcode sequences to the reference file
        ordering = seqs.Count[,1] %in% reffileseq
        # Count up how many sequences matching known barcodes were found in this primer
        # And add this number to the sequence stats table
        sequencestats[i,4] = sum(ordering * seqs.Count[,2])
        # Calculate the proportions of each barcode in this primer
        seqs.Count[,3] = seqs.Count[,2]/(ordering * sum(seqs.Count[,2]))
        # Make another column to add barcode names
        seqs.Count[,4] = ""
        # Match barcode sequences to the reference file and add them to the fourth column
        seqs.Count[ordering,4] = reffilename[match(seqs.Count[ordering,1], reffileseq)]
        # Reorder the columns
        seqs.Count = seqs.Count[,c(4, 1, 2, 3)]
        # And sort the barcodes by their proportions in the sample
        seqs.Count = seqs.Count[order(-ordering, -seqs.Count[,3]),]
        
        # If sum(ordering) is less than nrow(seqs.Count)
        # this means not every barcode in seqs.Count matched to the reference list
        # and, therefore there are "unique" barcodes for us to take into account/assign a unique label to
        if(sum(ordering) < nrow(seqs.Count)){
          
          # We want our proportion calculations to only include known barcodes, so set the unique's proportions to zero
          seqs.Count[((sum(ordering) + 1):(nrow(seqs.Count))),4] = 0
          # Assign each unique a name
          seqs.Count[((sum(ordering) + 1):(length(seqs.Count[,1]))),1] = paste("Unique", 1:sum(!ordering))
          
        }
        
        # If there are known barcodes in this sample
        if (sum(ordering) > 0){
          
          # Calculate the proportions of each known barcode
          seqs.Count[1:sum(ordering),4] = seqs.Count[1:sum(ordering),3]/sum(seqs.Count[1:sum(ordering),3])
          
        }
        
        # We now have a worked-up table of sequences and counts
        # First column is the barcode name, second is the actual sequence, third is the raw counts,
        # and fourth is the proportion in the sample
        colnames(seqs.Count)=c("Barcode","Sequence","Counts","Proportion")
        
        # Print to console what primer we found these reads in
        #print(paste("Barcodes found in primer", i))
        
        # Add these reads to our count of sequences found for this primer
        sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
        
        # If no reads found for this primer, don't bother adding to our final matrix
        if (nrow(seqs.Count) < 2){
  
          print("Didn't find any reads in this sample (checkpoint 5)")
          print(i)
  
          # But still count how many sequences were found in this sample
          # sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
  
          # next
        }
        if (length(seqs.Count) < 2){
  
          print("Didn't find any reads in this sample (checkpoint 6)")
          print(i)
  
          # But still count how many sequences were found in this sample
          sequencestats[i,2] = sum(sum(index), sequencestats[i,2])
  
          # next
        }
        
        # Overwrite this primer's data in the big table with the updated count list
        barcodeunique.hold[6:(nrow(seqs.Count) + 5),(8*i-7):(8*i-4)] = seqs.Count
        # And label the sample with the primer number
        barcodeunique.hold[1,(8*i-7):(8*i-2)] = rep(primer.list.cleaned[i,1], 6)
      }
      
      # Keep track of how many sequences we've analyzed so far
      sum = sum + length(fq)
      # And update the progress bar in the browser
      #incProgress(amount = 0, message = paste(sum/10^6, "Million Sequences Analyzed"))
      print(paste(sum/10^6, "million sequences analyzed"))
      
    }
    
    #### Post-processing ############################################################################################
    unprocessed <- barcodeunique.hold
    # Update the progress bar once we've analyzed all sequences
    
    barcodeunique.hold <- unprocessed
    #### 1/8/20 updated postprocessing  with 6 categories ####
    
    # For every sample
    for(i in 1:nrow(primer.list.cleaned)){
      
      # Pull the current list of sequences and counts for this animal
      sample <- barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-4)]
      
      # Omit empty rows
      sample <- na.omit(sample)
      # Convert numbers to numeric
      sample[,3] <- as.numeric(sample[,3])
      sample[,4] <- as.numeric(sample[,4])
      # Remove singletons
      sample <- sample[which(!sample[,3] < 2),]
      
      # These barcodes that were extracted & counted at least twice
      # Are what we count for the total reads extracted per primer
      sequencestats[i,3] <- sum(sample[,3])
      # And recalculate barcodes to include everything that passed
      sample[,4] <- sample[,3]/sum(sample[,3])
      
      # Check again to make sure there are reads in this sample left to analyze
      if(is.na(sample[1,1])){next}
      
      # Find the uniques (barcodes that don't have names/weren't found during stock characterization)
      uniques <- which(grepl("Unique", sample[,1]))
      
      # Pull data for only the known/named barcodes
      compare <- sample[-(uniques),]
      
      # And add stats about how many known barcodes we found to the stats table
      sequencestats[i,4] <- sum(compare[,3])
      
      # Proportion uniques is 1 - (the percentage of extracted sequences that were named barcodes)
      # We include this proportion on the sample data sheet as a quality control indicator
      barcodeunique.hold[3, (8*i-4)] <- 1-(sequencestats[i,4]/sequencestats[i,3])
      
      #### QC on named barcodes ####
      
      # Make sure there are enough barcodes to make a comparison, and that user wants the hamming distance flag
      if ((nrow(compare) > 1) & (mindist > 0)){
        
        print(paste("Starting regular hamming distance checks for sample", i))
        
        # Clear hamming distance cutoff point
        if(exists("hamming_space")){
          rm(hamming_space)
        }
        
        # Starting with the lowest barcode, find the barcode it's closest to
        for(j in 1:(nrow(compare)-1)){
          # Index of the barcode most similar to this one
          # 6/16/20: Changed back to hamming distance, not pairwise alignment
          w <- which.min(stringdist(compare[1:(nrow(compare)-j),2], compare[(nrow(compare)-j+1),2], method = "lv"))
          # And the barcode name
          best.fit <- compare[w,1]
          # And their distance apart
          num.dif <- min(stringdist(compare[1:(nrow(compare)-j),2], compare[(nrow(compare)-j+1),2], method = "lv"))
          # And make a flag for this barcode
          compare[(nrow(compare)-j+1),5] <- paste0(num.dif, "-", best.fit)
        }
        
        # After calculating distances, need to clean stuff up
        # Only want to keep flags that are at our min dist or less
        dists <- as.numeric(gsub("\\-.*", "", compare[,5]))
        compare[which((dists > mindist)|is.na(dists)),5] <- "NA-NA"
        
        # Let's find the parents for each potential offspring
        parents <- unlist(lapply(gsub("\\d\\-", "", compare[,5]), function(x){
          if(TRUE %in% grepl(paste0(x, "$"), compare[,1])){
            which(grepl(paste0(x, "$"), compare[,1]))
          }else{
            return(NA)
          }}))
        
        # This is the ratio of parent to offspring
        if(TRUE %in% (compare[, 3]/compare[parents,3] > 0.05)){
          
          # If the "offspring" is more than 5% of the size of the "parent"
          # It's probably real and not PCR error, so we can remove the flags
          compare[which(compare[, 3]/compare[parents,3] > 0.05),5] <- "NA-NA"
          parents[which(compare[, 3]/compare[parents,3] > 0.05)] <- NA
         
        }
        if(TRUE %in% (compare[,3]/compare[parents,3] >= 0.01)){
          # If the offspring and parent are close in proportion
          # Make sure the "parent" doesn't already have a hamming dist flag
          if(TRUE %in% !grepl("\\d-", compare[parents[which(compare[,3]/compare[parents,3] >= 0.01)],5])){
            w <- which(!grepl("\\d-", compare[parents[which(compare[,3]/compare[parents,3] >= 0.01)],5]))
            w <- which(compare[,3]/compare[parents,3] >= 0.01)[w]
            
            # If we don't have an ambiguous category yet
            # Move these barcodes to the top of the flagged section
            # These are barcodes with no flag
            unflagged <- which(grepl("NA", compare[,5]))
            # These are barcodes with a hamming dist flag
            flagged <- which(!grepl("NA", compare[,5]))
            # These are the flagged barcodes we want at the top
            #compare[w,]
            # So let's separate them out from our normal flagged barcodes
            flagged <- flagged[!(flagged %in% w)]
            # And finally let's reorder everything with the ambiguous barcodes in between
            compare <- compare[c(unflagged, w, flagged),]
            # And add a space between the ambiguous and almost certainly error barcodes
            hamming_space <- length(w)+1
            # Here are the "ambiguous" flagged barcodes for future reference
            compare[which(!grepl("NA", compare[,5]))[1:hamming_space],]
            }
        }else{
          unflagged <- which(grepl("NA", compare[,5]))
          # These are barcodes with a hamming dist flag
          flagged <- which(!grepl("NA", compare[,5]))
          compare <- compare[c(unflagged, flagged),]
      }
        
        # Let's rearrange to have our hamming distance flagged barcodes lower down the list
        compare <- compare[order(grepl("NA", compare[,5]), decreasing = TRUE),]
        
        # And add this back to the full data for this sample
        sample[,5] <- NA
        sample <- rbind(compare, sample[uniques,])
        
      }else{
        # If user doesn't want hamming distance flag, make that row all NAs
        sample[,5] <- NA
      }
      
      # Now we want to determine our initial list of likely real named barcodes for this sample
      # We will add an empty line in the data to indicate barcodes below the line are not trusted
      space <- NA
      
      # Index of named barcodes that have no hamming distance flag
      unflagged <- which(!grepl("\\d", sample[,5]) & !grepl("Unique", sample[,1]))
      
      # Now, apply input cutoff if selected (delete low uniques and find point where we stop trusting barcodes)
      if(cutoff == "Yes"){
        
        # If input is less than 200 then we default to a 1/200 cutoff
        # Otherwise cutoff is 1/input
        if(primer.list.cleaned[i,4] == Inf){
          cutoff_prop <- 1/10000
        }else{
          cutoff_prop <- min(1/as.numeric(primer.list.cleaned[i,4]), 1/200)
        }
        
        # 4/3/20: Removing this section for Taina's burst size analysis
        # # Delete uniques that are at a proportion less than cutoff (based on total # reads in our list)
        # if(nrow(sample[-intersect(which(sample[,3]/sum(sample[,3]) < cutoff_prop), uniques),]) > 0){
        #   
        #   sample <- sample[-intersect(which(sample[,3]/sum(sample[,3]) < cutoff_prop), uniques),]
        # }
        
        # Find the point at which the unflagged barcodes pass below 1/10,0000
        space <- which(sample[unflagged,3]/sum(sample[,3]) < cutoff_prop)[1]
        # Or the bottom of the trusted barcodes if everything passes cutoff
        if(is.na(space) & (nrow(sample) != length(uniques))){
          space <- max(unflagged)+1
          }
        
      }
      
      #### QC on uniques ####
      
      # Let's locate all of our uniques again
      uniques <- which(grepl("Unique", sample[,1]))
      # Here we will apply a cutoff based on the lowest proportion trusted named barcode
      if((TRUE %in% (sample[uniques,3] > sample[space-1,3]*100)) | (TRUE %in% (sample[uniques,4] > 0.01))){
        
        ## This is where we check for across stock contamination.
        trustedUniques <- sample[uniques[which((sample[uniques,3] > sample[space-1,3]*100) | (sample[uniques,4] > 0.01))],]
        hamDistUniques <- NA
        if (ContamCheck == T){
          WrongStockMatch <- which(trustedUniques[,2] %in% combinedSeqs)
          if(length(WrongStockMatch) > 0){
            # filter out any overlap between stocks by removing anything in M2set
            WrongStockBarcodes <- trustedUniques[WrongStockMatch,2]
            trustedUniques <- trustedUniques[-WrongStockMatch,]
          }
          if(nrow(trustedUniques) > 0){
            for(row in 1:nrow(trustedUniques)){
              # Index of the barcode most similar to this one
              # 6/16/20: Changed back to hamming distance, not pairwise alignment
              w <- which.min(stringdist(compare[,2], trustedUniques[row,2], method = "lv"))
              # And the barcode name
              best.fit <- compare[w,1]
              # And their distance apart
              num.dif <- adist(compare[w,2], trustedUniques[row,2])
              # And make a flag for this barcode
              trustedUniques[row,5] <- paste0(num.dif, "-", best.fit)
            }
            hammDist <- as.numeric(gsub("-.*","", trustedUniques[,5]))
            trustedUniques[which(hammDist > 3),5] <- "NA-NA"
            crossStockHamCheck <- sapply(trustedUniques[which(trustedUniques[,5] == "NA-NA"),2], function(x){
              SequenceOfInterest <- x
              minAlignStock <- adist(reffileseq, SequenceOfInterest)
              minAlignOthers <- adist(combinedSeqs, SequenceOfInterest)
              if(min(minAlignStock) < min(minAlignOthers)){
                matches <- reffilename[which(minAlignStock == min(minAlignStock))]
                matches <- paste0(adist(SequenceOfInterest, reffileseq[which.min(minAlignStock)]), "-", matches)
              }else{
                matches <- combinedNames[which(minAlignOthers == min(minAlignOthers))]
                matches <- paste0(adist(SequenceOfInterest, combinedSeqs[which.min(minAlignOthers)]), "-", matches)
              }
              return(paste(matches, collapse =","))
            })
            hammDist <- as.numeric(gsub("-.*","", crossStockHamCheck))
            crossStockHamCheck[which(hammDist > 3)] <- "NA-NA"
            trustedUniques[which(trustedUniques[,5] == "NA-NA"),5] <- crossStockHamCheck
            hamIndex <- grep(paste(incorrectStocks, collapse = "|"), trustedUniques[,5])
            if(length(hamIndex) > 0){
              hamDistUniques <- trustedUniques[hamIndex,]
              trustedUniques <- trustedUniques[-hamIndex,]
            }
          }


        }
        # If there are uniques just as high as trusted barcodes, move them up to the real
        # Pull our original list of named barcodes
        compare <- sample[!grepl("Unique", sample[,1]),]
        # And then add these uniques
        if(!is.na(space) & (space > 1) & nrow(trustedUniques)> 0){
          # If there are real trusted barcodes
          # Pull the uniques that are just as high as these real barcodes
          compare <- rbind(compare, trustedUniques)
        }
        # else{
        #   # If there are no trusted barcodes yet, just pull all the uniques
        #   # compare <- rbind(compare, sample[uniques,])
        # }
        
        
        # Reorder everything
        compare <- compare[order(compare[,3], decreasing = TRUE),]
        
        # 1/30/20: changing hamming distance checks for uniques specifically
        # To account for uniques that are slightly closer to another unique than to the main barcode
        # (Allowing it to pass the 5% cutoff)
        
        # Make sure there are enough barcodes to make a comparison, and that user wants the hamming distance flag
        if ((nrow(compare) > 1) & (mindist > 0)){
          
          print(paste("Starting unique hamming distance checks for sample", i))
          
          # Give the very first barcode an NA-NA since it's the highest
          compare[1, 5] <- "NA-NA"
          
          # Starting with the second highest barcode, find the lower barcode it's closest to
          for(j in 1:(nrow(compare)-1)){
            
            # Index of the barcode most similar to this one
            # 6/16/20: Changed back to hamming distance, not pairwise alignment
            w <- which.min(stringdist(compare[1:j,2], compare[(j+1),2], method = "lv"))
            # And the barcode name
            best.fit <- compare[w,1]
            # And their distance apart
            num.dif <- min(stringdist(compare[1:j,2], compare[(j+1),2], method = "lv"))
            # If these barcodes are more distant that our chosen distance, then don't flag
            if(num.dif > mindist){
              compare[j+1, 5] <- "NA-NA"
              next
            }
            # And make a flag for this barcode
            compare[(j+1),5] <- paste0(num.dif, "-", best.fit)
            
            # Next up, need to check our cutoffs
            if(compare[j+1, 3]/compare[w, 3] > 0.05){
              # If this barcode is more than 5% of its parent
              # Then it's probably not a PCR error, and therefore actually its own barcode
              compare[j+1, 5] <- "NA-NA"
            }
            
            #print(compare[j + 1,])
            compare[j+1, 3]/compare[w, 3]
            
          }
          # After calculating distances, need to clean stuff up
          
          # Let's find the parents for each potential offspring
          parents <- unlist(lapply(gsub("\\d\\-", "", compare[,5]), function(x){
            if(TRUE %in% grepl(paste0(x, "$"), compare[,1])){
              which(grepl(paste0(x, "$"), compare[,1]))
            }else{
              return(NA)
            }}))
          
          # This is the ratio of parent to offspring
          if(TRUE %in% (compare[,3]/compare[parents,3] >= 0.01)){
            # If there are flagged barcodes of similar proportions
            # Then first make sure they're above cutoff and need to be compared
            if(TRUE %in% (compare[which(compare[,3]/compare[parents,3] >= 0.01),4] > cutoff_prop)){
              
              # Identify the reads that meet this criteria
              w <- which(compare[which(compare[,3]/compare[parents,3] >= 0.01),4] > cutoff_prop)
              w <- which(compare[,3]/compare[parents,3] >= 0.01)[w]
              #compare[w,]
              
              # Next, see whether "parents" are offspring of another barcode
              # If there are any that aren't, we can take a look at them
              if(TRUE %in% !grepl("\\d-", compare[parents[w],5])){
                # Here are the barcodes that pass the cutoff
                # And have parents with no hamming distance flag
                w <- w[which(!grepl("\\d-", compare[parents[w],5]))]
                #compare[w,]
                
                # If we don't have an ambiguous category yet
                # Move these barcodes to the top of the flagged section
                unflagged <- which(grepl("NA", compare[,5]))
                #compare[w,]
                flagged <- which(!grepl("NA", compare[,5]))
                flagged <- flagged[!(flagged %in% w)]
                compare <- compare[c(unflagged, w, flagged),]
                # And add a space between the ambiguous and almost certainly error barcodes
                hamming_space <- length(w)+1
                # Here are the "ambiguous" flagged barcodes for future reference
                #compare[which(!grepl("NA", compare[,5]))[1:hamming_space],]
                
              }
              
            }else{
              # If none of these barcodes pass our cutoff, then they aren't ambiguous
              # They're definitely just not real
              if(exists("hamming_space")){rm(hamming_space)}
              if(FALSE %in% (compare[,3]/compare[parents,3] >= 0.01)){
                removeCutoff <- which(compare[,3]/compare[parents,3] < 0.01)
                compare <- compare[-removeCutoff,]
              }
            }
          }else{
            if(FALSE %in% (compare[,3]/compare[parents,3] >= 0.01)){
              removeCutoff <- which(compare[,3]/compare[parents,3] < 0.01)
              compare <- compare[-removeCutoff,]
            }
          }
          
          # If these samples aren't already separated into ambiguous and erroneous
          # Then just sort to have our hamming distance flagged barcodes lower down the list
          if(!exists("hamming_space")){
            compare <- compare[order(grepl("NA", compare[,5]), decreasing = TRUE),]
          }
          uniqueDistRetrieval <- which(compare[,2] %in% trustedUniques[,2])
          for(uniqueIndex in uniqueDistRetrieval){
            hamIndex <- which(trustedUniques[,2] == compare[uniqueIndex,2])
            compare[uniqueIndex,5] <- trustedUniques[hamIndex,5]
          }
          uniqueIndex <- intersect(grep("Unique", compare[,1]), which(compare[,5] == "NA-NA"))
          uniqueHams <- sapply(uniqueIndex, function(index){
            seq <- compare[index,2]
            minAlignStock <- adist(reffileseq, seq)
            matches <- reffilename[which(minAlignStock == min(minAlignStock))]
            matches <- paste0(adist(seq, reffileseq[which(minAlignStock == min(minAlignStock))]), "-", matches, collapse = ",")
          })
          compare[uniqueIndex,5] <- uniqueHams
          # And add this back to the full data for this sample
          sample <- rbind(compare, sample[which(!(sample[,1] %in% compare[,1])),])
          
          # And update our space indicating cutoff point
          space <- (which(compare[,4] < cutoff_prop)[1])
          # If nothing unflagged below cutoff, then just find point between flagged and unflagged
          if(is.na(space)){space <- which(grepl("\\d-", compare[,5]))[1]}
          # And if nothing flagged, just find the end of the "real" barcodes
          if(is.na(space)){space <- nrow(compare)+1}
          if (ContamCheck == T){
            WrongStockMatch <- which(sample[,2] %in% combinedSeqs)
            # filter out any overlap between stocks by removing anything in combinedStocks
            matchedBars <- sapply(sample[WrongStockMatch,2], function(z){
              paste(combinedNames[which(combinedSeqs == z)], collapse = ",")
            })
            sample[WrongStockMatch,5] <- matchedBars
            if(!is.null(nrow(hamDistUniques))){
              for(z in 1:length(hamDistUniques[,2])){
                sampleIndex <- which(sample[,2] == hamDistUniques[z,2])
                sample[sampleIndex,5] <- hamDistUniques[z,5]
              }
            }
          }
    
        }else{
          # If user doesn't want hamming distance flag, make that row all NAs
          sample[,5] <- NA
        }
        
      }
      
      #### Sample output formatting ####
      
      # Update the number of sequences matching known barcode to only include things passing our cutoff and QC
      if(is.na(space)){
        
        # If everything passes cutoff, just update the count to not include singletons we filtered out
        sequencestats[i,4] <- sum(sample[,3])
        
      } else if(space > 1){
        
        # Count only reads that passed our cutoff
        sequencestats[i,4] <- sum(sample[1:(space-1),3])
        
      } else if(space == 1){
        
        # If nothing passed our cutoff, set it to zero
        sequencestats[i,4] <- 0
        
      }
      
      # Lastly, let's separate out our different groups of data
      if(!is.na(space)){
        if(nrow(sample) >= space){
          
          # Everything that passed all the QC checks is the first group, that's the space
          if(space != 1){
            unflagged <- 1:(space-1)
          # Let's recalculate the proportions for these totally trusted barcodes
          sample[,4] <- c(sample[1:(space-1),3]/sum(sample[1:(space-1),3]), sample[space:nrow(sample),4])
          
          # Next let's identify our hamming distance flagged barcodes
          flagged <- which(grepl("\\d-", sample[, 5]))
          flagged <- flagged[-which(flagged %in% grep("Unique", sample[,1]))]
          # And let's move our unflagged named barcodes below the cutoff down with the uniques
          uniques <- c(1:nrow(sample))
          uniques <- uniques[-c(unflagged, flagged)]
          uniques <- uniques[order(sample[uniques, 3], decreasing = TRUE)]
          
          # And finally let's put this all back together in order
          sample <- sample[c(unflagged, flagged, uniques),]
          # And add our space to show the point where we stop trusting barcodes
          sample <- rbind(sample[1:(space-1),], NA, sample[space:nrow(sample),])
          
          # If there's an "ambiguous" section then we add a space to indicate that too
          if(exists("hamming_space")){
            unflagged <- 1:(space-1)
            ambiguous <- space:(space+hamming_space-1)
            flagged <- (space+hamming_space):nrow(sample)
            sample <- rbind(sample[unflagged,], sample[ambiguous,],NA, sample[flagged,])
          }
        }else{
            sample <- rbind(NA, sample)
        }
      }
    }
      # (Don't need to rearrange if everything from the start passed QC)
      # And now we can write over the old sample data for this primer with our QC'd version
      barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-3)] <- NA
      barcodeunique.hold[6:(nrow(sample) + 5), (8*i-7):(8*i-3)] <- sample
      ## add a percent contamination (percent of reads marked as incorrect stock)
      if(ContamCheck == TRUE){
      hamIndex <- which(sample[,2] %in% combinedSeqs)
      propContam <- sum(as.numeric(sample[hamIndex,3]))/sum(as.numeric(sample[!is.na(sample[,3]),3]))
      barcodeunique.hold[3, (8*i-3)] <- propContam
      }
    }
    
    # Once we've applied QC to each individual sample
    # We can check for index hopping between samples
    
    # Eventually I want to change the app to keep track of which samples are from which animals
    # So that we can have runs with multiple samples from multiple different animals
    # Without worrying about stuff getting incorrectly flagged when it's in two samples from the same animal
    # (This is probably a rare occurrence anyway though)
    
    pre_contam <- barcodeunique.hold
    barcodeunique.hold <- pre_contam
    
    #### Index hopping checker ####
    # With multiplexed runs we expect a small amount (0.1-0.5%) of reads to bleed over from other samples
    # This is called index hopping or index misassignment
    # For this reason, we want to compare barcodes between our multiplexed samples
    # And if a barcode is very high in one sample and very low in another
    # These might be index-hopped reads and not really in the lower sample
    
    # Build a matrix of all data in this run
    # With counts for all barcodes that we considered significant
    
    # Then we can make our final list of trusted barcodes and calculate final proportions and sample stats from them
    
    # If user set an index-hopping threshold
    if(contam < 100){
      
      # This will be our list of real barcodes to check
      if(!is.na(reffilepath)){
        # Make a table of all barcodes in the reference file and their names
        check <- data.frame(name = ShortRead::id(reffile), sequence = sread(reffile))
      }else{
        check <- data.frame(name = NA, sequence = NA)
      }
      
      # This will be our list of all reads to compare against
      all <- check
      
      # Update 1/31/20: want to check only "real" barcodes for index-hopping as always
      # But want to compare those real barcodes against all barcodes, real or error
      
      # Pull data for each sample
      for(i in 1:nrow(primer.list.cleaned)){
        
        # Pull the current list of sequences and counts for this animal
        sample <- barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-6):(8*i-5)]
        
        # Update formatting
        colnames(sample) <- c("sequence", primer.list.cleaned[i,1])
        sample[,2] <- as.numeric(sample[,2])
        
        # For full list, pull all data from this sample
        all <- merge(all, na.omit(sample), by = "sequence", all = TRUE)
        
        # Find the cutoff point
        space <- which(is.na(sample))[1]
        
        if(is.na(space)){
          
          print("no cutoff!")
          print(i)
          
        } else{
          
          # And remove everything below cutoff
          sample <- sample[1:(which(is.na(sample))[1]-1),]
        }
        
        # Then add this sample to our table of all data in the run
        check <- merge(check, sample, by = "sequence", all = TRUE)
        
      }
      
      # Remove rows that weren't counted at all
      check <- check[rowSums(check[,3:ncol(check)], na.rm = TRUE) > 0,]
      # The full list has an order of magnitude more barcodes than the real
      # So we can also take all those missing barcodes out of the full list
      all <- all[all[,1] %in% check[,1],]
      # And make sure there's no duplicates
      all <- unique(all)
      
      # Check whether this run has animal names
      if(length(unique(runinfo$Animal[!is.na(runinfo$Animal)])) > 1){
        
        # Need to associate primers and animal names
        # All of the primers
        #colnames(all[,3:ncol(all)])
        # All of the animal names & associated primers
        #data.frame(runinfo$Animal[!is.na(runinfo$Animal)], runinfo$Barcodes[!is.na(runinfo$Animal)])
      
        # Check to make sure the primers are listed in numerical order in our runinfo file
        # Set primer_prefix based on fwdrev
        if (fwdrev == "Forward") {
          primer_prefix = "INT."
        } else if (fwdrev == "Reversed") {
          primer_prefix = "VPX."
        } else {
          stop("fwdrev not set")
        }
        runinfo$Barcodes[!is.na(runinfo$Animal)] <- toupper(runinfo$Barcodes[!is.na(runinfo$Animal)])
        if(FALSE %in% (colnames(all[,3:ncol(all)]) == gsub(toupper(primer_prefix), "", runinfo$Barcodes[!is.na(runinfo$Animal)]))){
          stop("Implement formatting for primers listed out of order, or atypical primer prefix used in runinfo file (VPX or Int?)")
        }else{
          #colnames(all)[3:ncol(all)] <- runinfo$Animal[!is.na(runinfo$Animal)]
          colnames(check)[3:ncol(check)] <- runinfo$Animal[!is.na(runinfo$Animal)]
        }
        
        hold <- check
        
        # Change the primers to animal names
        
        # For each barcode, find the samples that came from a different sample
        for(i in 3:ncol(check)){
          
          # Get the name of the animal for this sample
          x <- gsub("\\..+", "", colnames(check)[i])
          # Sum up the counts for that barcode in samples from other animals
          #x <- colnames(check)[4]
          #colnames(all)[which(!grepl(x, colnames(all)[3:ncol(all)]))+2]
          
          # And calculate the ratio of the counts to the sum of counts for that barcode in other animals
          hold[,i] <- check[,i]/apply(all[, which(!grepl(x, colnames(check)[3:ncol(check)]))+2, drop = F], 1, function(y){
            # If there are no counts for this barcode from another animal, return one
            if(sum(!is.na(y))==0){return(1)}
            # Otherwise return the max counts
            else{return(sum(y, na.rm = TRUE))}
          })
        }
        
        check <- hold
        
      }else{
        
        # For each barcode, find the ratio of the counts to the maximum count
        check[,3:ncol(check)] <- check[,3:ncol(check)]/apply(all[, 3:ncol(all)], 1, sum, na.rm = TRUE)
        
      }
      
      if(nrow(check) > 0){
        # And find where this ratio is less that 0.002
        hopped <- which(check < 0.002, arr.ind = TRUE)
      }
      
      for(i in 1:nrow(primer.list.cleaned)){
        
        # Pull the current list of sequences and counts for this animal
        sample <- barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-2)]
        
        if(nrow(check) > 0){
          # And flag barcodes that we think might have index hopped from a different sample
          sample[which(sample[,2] %in% check[hopped[(hopped[,2] == (i + 2)),1],1]),6] <- "YES"
        }
        
        # Find the cutoff point
        space <- which(is.na(sample))[1]
        
        # Reorder to have the newly flagged barcodes further down the list
        if(is.na(space)){
          
          print("no cutoff!")
          print(i)
          
        } else if(space == 1){
          # If everything is below cutoff, don't change anything
          next
          
        } else{
          
          # Pull the list of trusted barcodes
          unflagged <- sample[(1:(space-1)),]
          # Then take the newly flagged barcodes and add them to the flagged barcode list
          flagged <- rbind(unflagged[!is.na(unflagged[,6]),], sample[(space+1):nrow(sample),])
          # And remove these flagged barcodes from our trusted barcode list
          unflagged <- unflagged[is.na(unflagged[,6]),]
          
          # Recalculate proportions to only include our trusted barcodes
          unflagged[,3] <- as.numeric(unflagged[,3])
          unflagged[,4] <- unflagged[,3]/sum(unflagged[,3])
          
          # And mark barcodes that passed our index hopping check
          if(nrow(unflagged) > 0){unflagged[,6] <- "NO"}
          
          # Finally, reorder the sample to have the flagged barcodes lower down the list
          sample <- rbind(unflagged, NA, flagged)
          sample2 <- sample[-which(duplicated(sample)),]
        }
        
        # If there's an ambiguous category for this data,
        # move the index hopped barcodes below this too
        
        # We can determine this by figuring out if the second and third empty lines are distant
        if(((which(is.na(sample))[3]-which(is.na(sample))[2]) > 1) & ("YES" %in% sample[,6])){
          
          # Find where the second space is
          hamming_space <- which(is.na(sample))[2]
          
          # And move the index hopped samples below this space
          
          # First pull everything above the second space
          unflagged <- 1:(hamming_space-1)
          # And remove the index hopped barcodes
          #which(sample[,6] == "YES")
          unflagged <- unflagged[!(unflagged %in% which(sample[,6] == "YES"))]
          # Then pull everything else
          flagged <- (hamming_space+1):nrow(sample)
          # And put the hopped barcodes in between
          sample <- rbind(sample[unflagged,], NA, sample[which(sample[,6] == "YES"),], sample[flagged,])
          
        }
        
        # Update the sample in the main table
        barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-2)] <- sample
        
      }
      
    }
    
    #### Short barcode flag ####
    # 6/17/20: Adding flag for short barcodes
    # At least for now not going to reorder anything, just ID the short ones
    # They can still be real barcodes even though they're short
    # Deletion should happen naturally over time
    
    print("Checking for short barcodes")
    for(i in 1:nrow(primer.list.cleaned)){
      
      # Pull the data for this sample
      sample <- barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-1)]
      # Let's just look at the barcode sequences
      barcodes <- na.omit(sample[,2])
      # Check against several known short barcodes, with the deletion in different places
      # Short barcode 1
      w <- which(pairwiseAlignment(barcodes, "CCCCTCCCCCTCCAGGACTAGCATAAACGCGCGT", scoreOnly = T) > 0)
      sample[which(sample[,2] %in% barcodes[w]),7] <- "short"
      # Short barcode 2
      w <- which(pairwiseAlignment(barcodes, "GACCTCCTCCTCCTCCCCCTCCAGGACTAGCAGT", scoreOnly = T) > 0)
      sample[which(sample[,2] %in% barcodes[w]),7] <- "short"
      # Short barcode 3
      w <- which(pairwiseAlignment(barcodes, "GAGACCAGGACCTCCTCCTCCTCCCCCTCCAGGA", scoreOnly = T) > 0)
      sample[which(sample[,2] %in% barcodes[w]),7] <- "short"
      # Short barcode 4
      w <- which(pairwiseAlignment(barcodes, "CCCTCCTCCCCCTCCAGGACTAGCATAAACGCGT", scoreOnly = T) > 0)
      sample[which(sample[,2] %in% barcodes[w]),7] <- "short"
      
      # Now add this to our output data file
      barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-1)] <- sample
      
    }
    
    # Formatting after QC
    #Now that we've flagged barcodes, checked our uniques, applied our cutoff, and check for index hopping
    # Let's update the summary statistics for each sample
    # And remove samples from our output file that didn't find any reads
    
    for(i in 1:nrow(primer.list.cleaned)){
      
      # Pull the data for this sample
      sample <- barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-1)]
      # And find the cutoff point
      
      # Find the cutoff point
      space <- which(is.na(sample))[1]
      
      # Reorder to have the newly flagged barcodes further down the list
      if(is.na(space)){
        
        print("no cutoff!")
        print(i)
        
      } else if(space == 1){
        # If everything is below cutoff, trusted barcode count is zero
        sequencestats[i,4:5] <- 0
        
      }else{
        # Recalculate read and barcode counts to only include trusted barcodes
        sequencestats[i,4] <- sum(as.numeric(sample[1:(space-1),3]))
        sequencestats[i,5] <- space-1
      }
      
      # 3/12/20
      # Adding hamming distance flags for real barcodes as well
      # Let's just pull any barcodes in our "real" category and check to see if they're related to anything higher
      compare <- sample[1:space-1,]
      
      # Make sure there are enough barcodes to make a comparison, and that user wants the hamming distance flag
      if ((nrow(compare) > 1) & (mindist > 0)){
        
        print(paste("Starting final hamming distance checks for sample", i))
        
        # Give the very first barcode an NA-NA since it's the highest
        compare[1, 5] <- "NA-NA"
        
        # Starting with the second highest barcode, find the lower barcode it's closest to
        for(j in 1:(nrow(compare)-1)){
          
          # First make sure this barcode isn't already matched up
          if(grepl("NA-NA", compare[j+1,5])){
            
            # Index of the barcode most similar to this one
            # 6/16/20: Changed back to hamming distance, not pairwise alignment
            w <- which.min(stringdist(compare[1:j,2], compare[(j+1),2], method = "lv"))
            # And the barcode name
            best.fit <- compare[w,1]
            # And their distance apart
            num.dif <- min(stringdist(compare[1:j,2], compare[(j+1),2], method = "lv"))
            # And make a flag for this barcode
            compare[(j+1),5] <- paste0(num.dif, "-", best.fit)
            # If these barcodes are more distant that our chosen distance, then don't flag
            if(num.dif > mindist){
              compare[(j+1), 5] <- "NA-NA"
              next
            }
            
            #stop("Found one!")
            # And make a flag for this barcode
            compare[(j+1),5] <- paste0(num.dif, "-", best.fit)
            
            
          }
          
        }
  
        # And add this new info back to the full data for this sample
        sample[1:space-1,] <- compare
        # Update the sample in the main table
        barcodeunique.hold[6:nrow(barcodeunique.hold), (8*i-7):(8*i-1)] <- sample
      }
      
    }
    
    # Add input values to stats table
    sequencestats[,6] <- primer.list.cleaned[,4]
    # Then add all of the sample stats to the summary at the top of each sample table
    barcodeunique.hold[3, seq(1, ncol(barcodeunique.hold), 8)] <- primer.list.cleaned[,4]
    barcodeunique.hold[3, seq(2, ncol(barcodeunique.hold), 8)] <- sequencestats[,4]
    barcodeunique.hold[3, seq(3, ncol(barcodeunique.hold), 8)] <- sequencestats[,5]
    
    #### Add to batch summary ####
    # Now that we've done all of the analysis for this run
    # Let's make some changes to account for the fact that we're doing multiple different runs
    
    # First I want to integrate the sequence stats into the runinfo like I usually do manually
    
    # Reorder sequencestats to match runinfo
    sequencestats <- sequencestats[order(unlist(lapply(sequencestats$Primer, function(x){
      grep(paste0(x, "$"), runinfo$Barcodes)
    }))),]
    # And add columns to runinfo
    
    # Abbey's version (deprecated w/ most recent R)
    ## runinfo[,(ncol(runinfo)+1):(ncol(runinfo)+ncol(sequencestats))] <- NA
    ## runinfo[1:nrow(sequencestats),(ncol(runinfo)-ncol(sequencestats)+1):(ncol(runinfo))] <- sequencestats
    ## colnames(runinfo)[(ncol(runinfo)-ncol(sequencestats)+1):ncol(runinfo)] <- colnames(sequencestats)
    
    # New version, don't hold me to this?
    tmp = merge(runinfo, sequencestats, by.x = 0, by.y = 0, all = T)
    tmp$Row.names = as.numeric(tmp$Row.names)
    tmp = tmp[order(tmp$Row.names),]
    tmp$CountsToInputRatio <- tmp$`# of Sequences Matching Known Barcode`/tmp$`Sequencing Input`
    tmp$Row.names = NULL
    rownames(tmp) = NULL
    
    runinfo = tmp
    tmp = NULL
    
    # And save runinfo to our list of excel sheets
    list_of_datasets <- append(list_of_datasets, list(runinfo))
    names(list_of_datasets)[length(list_of_datasets)] <- batchinfo$runinfo[file]
    
    #### Add sample names and dates to our summary tables ####
    
    runinfo$Date <- as.character(runinfo$Date)
    
    for(x in (primer.list.cleaned$`5' Index Name`)){
      barcodeunique.hold[1, grep(paste0(x, "$"), barcodeunique.hold[1,])] <- runinfo[grep(paste0(x, "$"), runinfo$Barcodes), 2:4]
    }
    
    # And remove samples that had no reads whatsoever from our main table
    # Their counts will still be included in the summary stats table
    i <- which(sequencestats[,3] == 0)
    barcodeunique.hold[,sort(outer(X = 8*i-7, Y = (0:7), FUN = "+"))] <- NULL
    
    #### And add the data to our all samples sheet ####
    if(!exists("all_samples")){
      all_samples <- barcodeunique.hold
    }else{
      # Make sure two dataframe are same length
      if(nrow(barcodeunique.hold) > nrow(all_samples)){
        all_samples[(nrow(all_samples)+1):nrow(barcodeunique.hold),] <- NA
      }else if(nrow(barcodeunique.hold) < nrow(all_samples)){
        barcodeunique.hold[(nrow(barcodeunique.hold)+1):nrow(all_samples),] <- NA
      }
      # Then add the new data
      all_samples[,(ncol(all_samples)+1):(ncol(all_samples)+ncol(barcodeunique.hold))] <- barcodeunique.hold
    }
  }
  
  #### Build Proportions Matrix (updated 12/19/19) ####
  
  # The main steps here are:
  # 1. Set up & add full data sheet
  # 2. Make & add the actual proportions matrix
  # 3. (iterate 1 & 2 as many times as need)
  # 4. Add batchinfo
  # 5. Name and write output file
  
  # If we are doing both the full and the normal then we need to output twice,
  # otherwise only need to do this once.
  if(full){
    noutputs <- 2
  }else{
    noutputs <- 1
  }
  
  # For each output file we're making
  for(output in 1:noutputs){
    
    # First pick which data we need to use
    if(output == 1){
      # First output is the normal one
      data_holder <- all_samples
      data_output <- list_of_datasets
    }else if(output == 2){
      # Second output is the "full" one (includes uniques in matrix)
      data_holder <- all_samples_full
      data_output <- list_of_datasets_full
    }
    
    # Then determine whether we need to order or separate data by animal
    
    # Get the list of animals in this run
    animal_list <- unique(as.character(data_holder[1, seq(1, ncol(data_holder), 8)]))
    
    # If there are more samples than animals, we want to loop through & format each animal
    if(length(animal_list) < ncol(data_holder)/8){
      
      for(animal in animal_list){
        
        # Pull all the data for this animal
        all_animal <- data_holder[,unlist(lapply(which(data_holder[1, seq(1, ncol(data_holder), 8)] == animal), function(x){((x*8)-7):((x*8))}))]
        # If there are multiple samples for this animal, then sort and build matrix
        if(ncol(all_animal)/8 > 1){
          # Sort data by sample date or time
          # If in hour format, order by hour
          if(TRUE %in% grepl("Hour", all_animal[1, seq(3, ncol(all_animal), 8)])){
            all_animal <- all_animal[, unlist(lapply(order(as.numeric(gsub("Hour ", "", all_animal[1, seq(3, ncol(all_animal), 8)]))), function(x){
              ((x*8)-7):((x*8)-1)
            }))]
          }else{
            # If in date format, order by date
            # First, check for and convert from different date formats
            all_animal[1, seq(3, ncol(all_animal), 8)] <- unlist(lapply(all_animal[1, seq(3, ncol(all_animal), 8)], function(x){
              if(suppressWarnings(!is.na(as.numeric(x)))){
                as.character(as.Date(as.numeric(x), origin = "1899-12-30"))
              }else if(!is.na(anydate(as.character(x)))){
                as.character(anydate(as.character(x)))
              }else if(!is.na(as.Date(as.character(x), tryFormats = "%m/%d/%y"))){
                as.character(as.Date(as.character(x), tryFormats = "%m/%d/%y"))
              }else if(is.na(x)){x}else{
                warning(paste0("Unable to interpret time/date: "), x)
                as.character(x)
              }
            }))
            
            # Want to also put plasma samples first
            # And DNA before cDNA
            # And also if multiple samples with same name, want to give a unique name
            headers <- all_animal[1, unlist(lapply(seq(3, ncol(all_animal), 8), function(x){(x-2):(x)}))]
            headers <- as.data.frame(t(matrix(headers, ncol = ncol(all_animal)/8)))
            headers[,4] <- 4
            
            # Assign priority for sample types
            headers[which(grepl("plasma", tolower(headers[,2]))), 4] <- 1
            headers[which(grepl("[^c]dna", tolower(headers[,2]))), 4] <- 2
            headers[which(grepl("cdna", tolower(headers[,2]))), 4] <- 3
            headers[which(grepl("rna", tolower(headers[,2]))), 4] <- 3
            headers[which(!grepl("plasma|[^c]dna|cdna|rna", tolower(headers[,2]))), 4] <- 4
            
            # Now we can sort by date first and sample second
            all_animal <- all_animal[, unlist(lapply(order(as.Date(unlist(headers[,3])), headers[,4]), function(x){
              ((x*8)-7):((x*8))
            }))]
            
            # Lastly, check and make sure all samples have unique names
            headers <- all_animal[1, unlist(lapply(seq(3, ncol(all_animal), 8), function(x){(x-2):(x)}))]
            headers <- as.data.frame(t(matrix(headers, ncol = ncol(all_animal)/8)))
            
            # If any names are duplicated, loop through and give each a unique name
            while(TRUE %in% duplicated(headers)){
              # Find which samples are duplicated
              dups <- which(unlist(headers[,3]) == unlist(headers[anyDuplicated(headers),3]))
              headers[dups, 2] <- paste(unique(headers[dups, 2]), 1:length(dups))
            }
            # And lastly add these names back to the datasheet
            all_animal[1, seq(2, ncol(all_animal), 8)] <- headers[,2]
          }
            
          # Remove NAs from headers to avoid Windows compatibility issues
          all_animal[1, is.na(all_animal[1,])] <- ""
          colnames(all_animal) <- all_animal[1,]
          all_animal <- all_animal[2:nrow(all_animal),]
          
          # Now that the animal datasheet is properly formatted, add to our output
          data_output <- append(data_output, list(all_animal))
          names(data_output)[length(data_output)] <- paste(animal, "samples")
          
          # And now that we have all the data saved, we can build a matrix
          
          # Use a reference list as backbone if we have one
          if(!is.na(reffilepath)){
            barcodes <- readFasta(reffilepath)
            barcodes <- data.frame(ShortRead::id(barcodes), sread(barcodes))
            colnames(barcodes) <- c("Name", "Sequence")
          }else{
            barcodes <- data.frame(Name = NA, Sequence = NA)
          }
          
          # Now we should be able to just paste in my usual proportion matrix code
          merge <- barcodes
          
          # So for each primer we want to pull only the sequences and counts *before* the space to indicate our cutoff
          # i = start + 2
          if(length(colnames(all_animal)) > 6){
            for(i in seq(3, length(colnames(all_animal)), 8)){
              # This should be our counts
              rep <- all_animal[,(i-1):i]
              # Find start of data and trim off headers
              rep <- rep[((which(is.na(rep[[1]]))[1])+2):length(rep[[1]]),]
              # Find the space that indicates cutoff and ignore things past that
              # If using readxl have to change this bit to use is.na because tidyverse/tibble stuff (ugh)
              rep <- rep[1:(which(is.na(rep[[1]]))[1]-1),]
              #print(length(rep[[1]]))
              
              rep[[2]] <- as.numeric(rep[[2]])
              rep <- na.omit(rep)
              colnames(rep) <- c("Sequence", paste(colnames(rep)[1], colnames(rep)[2]))
              
              merge <- merge(merge, rep, by = "Sequence", all = TRUE)
            }
            
            # Change all the NAs to zero
            merge[(is.na(merge))] <- 0
            
            # Remove barcodes that don't have any reads
            merge$sums <- rowSums(merge[,3:length(colnames(merge))])
            table <- merge[(merge$sums > 0), 1:(length(colnames(merge))-1)]
            
            # Calculate proportions
            for(i in 3:length(colnames(table))){
              table[[i]] <- table[[i]]/sum(table[[i]])
            }
            
            # Also nice to sort this by the first timepoint just to get in a meaningful order
            table <- table[(order(table[[3]], decreasing = TRUE)),]
            data_output <- append(data_output, list(table))
            names(data_output)[length(data_output)] <- paste(animal, "matrix")
          }
          
        }else{
          # If only one sample for this animal, just add that sample
          
          # Remove NAs from headers to avoid Windows compatibility issues
          all_animal[1, is.na(all_animal[1,])] <- ""
          colnames(all_animal) <- all_animal[1,]
          all_animal <- all_animal[2:nrow(all_animal),]
          
          # Now that the animal datasheet is properly formatted, add to our output
          data_output <- append(data_output, list(all_animal))
          names(data_output)[length(data_output)] <- paste(animal, "samples")
          
        }
        
      }
      
      
    }else{
      # If only one sample per animal, add just one big sheet to workbook
      
      # Remove NAs from headers to avoid Windows compatibility issues
      all_samples[1, is.na(all_samples[1,])] <- ""
      colnames(all_samples) <- all_samples[1,]
      all_samples <- all_samples[2:nrow(all_samples),]
      
      # Now that the datasheet is properly formatted, add to our output
      data_output <- append(data_output, list(all_samples))
      names(data_output)[length(data_output)] <- "All samples"
      
    }
    
    # Lastly, add settings for this run
    batchinfo[(nrow(batchinfo)+1):(nrow(batchinfo) + 12),] <- NA
    
    batchinfo[(nrow(batchinfo)-11):nrow(batchinfo), 1] <- c("App version", "Directory",
                                                            "Adapter sequence",
                                                            "Reference file path",
                                                            "Strandedness", "Reference sequence",
                                                            "Mismatches allowed",
                                                            "Location of barcode relative to reference",
                                                            "Length of barcode",
                                                            "Hamming distance to flag",
                                                            "Index hopping threshold",
                                                            "Analysis date")
    batchinfo[(nrow(batchinfo)-11):nrow(batchinfo), 2] <- c(version, folder, adapter, reffilepath,
                                                            fwdrev, ref, mismatches, direction,
                                                            numbases, mindist, contam, as.character(Sys.time()))
    
    data_output <- append(data_output, list(batchinfo))
    names(data_output)[length(data_output)] <- "Analysis settings"
    
    # And now we can write this to an excel workbook
    
    # If list names are more than the 31 allowed characters, truncate to preserve unique names
    names(data_output) <- substr(names(data_output), nchar(names(data_output))-30,
                                      nchar(names(data_output)))
    foldername <- unlist(str_split(folder, "/"))
    foldername <- foldername[length(foldername)-1]
    if(output == 1){
      name <- paste0(Sys.Date(), " ", foldername, " Analysis")
    }else if(output == 2){
      name <- paste0(Sys.Date(), " ", foldername, " Full Analysis")
    }
    
    # Write to excel workbook 
    wb <- createWorkbook()
    negStyle <- createStyle(fontColour = "#000000", bgFill = "#ffcccc")
    lapply(names(data_output), function(name){
      
      addWorksheet(wb, sheetName = name)
      writeData(wb, sheet = name, x = data_output[[name]])
    })
    setColWidths(wb, names(data_output)[1], cols = ncol(data_output[[1]]), widths = "auto")
    conditionalFormatting(wb, sheet= names(data_output)[1], cols = ncol(data_output[[1]]), rows = 2:nrow(data_output[[1]]),rule = "< 1", style = negStyle)
    if(file.exists(paste0(folder, "/", name, ".xlsx"))){
      
  
      # Use regex to see if multiple files have this name
      if(TRUE %in% grepl(paste0(name, " \\(\\d\\).xlsx"), list.files(folder))){
        
        # And if there are multiple, increment name by one
        saveWorkbook(wb, paste0(folder, "/", name, " (", sum(grepl(paste0(name, " \\(\\d\\).xlsx"), list.files(folder))) + 1, ").xlsx"), overwrite = F)
        
      }else{
        # If there's only one file with this name, just add a one to the name
        saveWorkbook(wb, paste0(folder, "/", name, " (1).xlsx"), overwrite = F)
      }
    }else{
      # Write like normal
      saveWorkbook(wb, paste0(folder, "/", name, ".xlsx"), overwrite = F)
    }
    
  }
}
