# This is the R part of the pipeline used to analyze NGS data received from Fasteris. Study
# species is brown bear, Ursus arctos.

library(parallel)
library(tidyr)

source("microsatTabExtract.R")

do.chunk.init <- FALSE
do.chunk1 <- FALSE
do.chunk2 <- FALSE
do.chunk3 <- FALSE
do.chunk4 <- TRUE
do.chunk5 <- TRUE

# parallel stuff
if (do.chunk.init) {
  # start cores
  if (Sys.info()["sysname"] == "Windows") {
    ncores <- 4
  } else {
    ncores <- 48
  }
  
  message(sprintf("Powering up %d processes.", ncores))
  cl <- makeCluster(ncores, outfile = "clusterfuck.txt")
  on.exit(stopCluster(cl))
}

# Store results of step 1 and 2 in this folder
dir.ngsfiler <- "./DAB/1_ngsfilters_hiseq1"
if (!dir.exists(dir.ngsfiler)) {
  dir.create(dir.ngsfiler)
}

dir.lsl <- "./DAB/3_lib_sample_locus_hiseq1" # notice no trailing slash
if (!dir.exists(dir.lsl)) { # create if doesn't exist
  dir.create(dir.lsl)
}

dir.uniq.tab <- "./DAB/2_uniq_tab_hiseq1"
if (!dir.exists(dir.uniq.tab)) {
  dir.create(dir.uniq.tab)
}

if (do.chunk1) {
  ####################################################
  #### Step 1: Create library-sample-locus files. ####
  ####################################################
  
  # Example:
  # R -q --vanilla --args MICROSAT.PCR_UA_MplexRout1_03.tab UARef < microsat_tabextract.r
  #                         |_ locus file with all samples    |_ sample name
  
  # to extract data for all loci of all samples, we need to create all
  # combinations
  message("Performing step 1.")
  
  # needed for this step
  clusterExport(cl = cl, varlist = c("microsatTabExtract"))
  clusterEvalQ(cl = cl, expr = {library(data.table)})
  
  # load data and sample names
  message("Importing files.")
  sn <- list.files(path = dir.ngsfilter, pattern = ".ngsfilter", full.names = TRUE)
  
  message(sprintf("Found %s ngs filters", length(sn)))
  names(sn) <- basename(sn)
  sn <- as.list(sn)
  sn <- sapply(sn, read.table, simplify = FALSE)
  
  inputfile <- list.files(dir.uniq.tab, pattern = "^MICROSAT.*\\.uniq\\.tab$", full.names = TRUE)
  libnum <- gsub("^.*_JFV-(\\d+)_UA_.*\\.uniq.tab$", "\\1", basename(inputfile))
  libnum <- sprintf("%02d", as.numeric(libnum))
  
  inputfile <- split(inputfile, f = libnum)
  
  # check that correct ngs filter (sample names) are assigned to correct library data
  message("Do input numbers match library ngsfilter file?")
  print(data.frame(ngsfilter = names(sn), input = names(inputfile)))
  
  # extract sample names
  samplenames <- sapply(sn, FUN = function(x) as.character(x$V2), simplify = FALSE)
  samplenames <- sapply(samplenames, FUN = function(x) {
    unique(gsub("^(.*)_[[:alnum:]]{3}_PP\\d+$", "\\1", x))
  }, simplify = FALSE)
  
  # Create a combination to extract alleles from all loci for a given sample.
  slc <- mapply(FUN = function(x.sample, x.input) {
    expand.grid(samplename = x.sample, inputfile = x.input)
  }, samplenames, inputfile, SIMPLIFY = FALSE)
  slc <- do.call(rbind, slc)
  slc$samplename <- as.character(slc$samplename)
  slc$inputfile <- as.character(slc$inputfile)
  rownames(slc) <- NULL
  
  message(sprintf("Processing %d files to %s.", nrow(slc), dir.lsl))
  
  # # Split the data into ncores chunks which are to be processed in parallel.  
  slc <- split(slc, f = slc$samplename)
  
  # progress report
  numsam <- length(slc)
  i <- 1
  
  out <- sapply(X = slc, FUN = function(x, dir.lsl) {
    message(sprintf("Processing %s (%d/%d)", unique(x$samplename), i, numsam)) # scoping out of current env!
    i <- i + 1
    out <- parApply(cl = cl, X = x, MARGIN = 1, FUN = function(y, outdir) {
      microsatTabExtract(filename = y["inputfile"], samplename = y["samplename"], outdir = outdir)
    }, outdir = outdir)
  }, outdir = dir.lsl)
  NULL
}

if (do.chunk2) {
  ######################################
  #### Step 2: Adds series to reads ####
  ######################################
  message("Performing step 2.")
  # load motif table (copied from de Barba et al. 2016).
  # This script assigns to each sequence a serie number that is added to the data in a new column.
  # It requires two arguments : 1) the tab filename (output from microsat_tabextract.r) 2) the sequence of the 
  # motif of the corresponding microsat marker
  
  # Example:
  # microsatTabToseries.py -f  MICROSAT.PCR_UARef_03.tab -m ctat
  #                                 |_ allele file            |_ motif from the paper
  
  # find files created by microsatTabExtract()
  # NOTICE: ZF locus is excluded from the analysis
  # NOTICE: copy microsatTabToseries.py to where the uniq.tab files are being processed
  # lsl = library-sample-locus files
  
  file.copy(from = "microsatTabToseries.py", to = paste(dir.lsl, "microsatTabToseries.py", sep = "/"),
            overwrite = TRUE)
  
  lsl <- list.files(path = dir.lsl, pattern = "MICROSAT.PCR_.*uniq.tab$")
  # load motif data
  motif <- read.table("locus_motifs.txt", header = TRUE,
                      colClasses = c("character", "character"))
  
  message(sprintf("Found %d files to process.", length(lsl)))
  print(head(lsl))
  
  if (!file.exists("microsatTabToseries.py")) stop("microsatTabToseries.py not found")
  
  # nodes for some reason are not set with current working dir, so do it explicitly
  # on each node
  clusterExport(cl = cl, varlist = "dir.lsl")
  clusterEvalQ(cl = cl, expr = {setwd(dir.lsl)})
  
  parSapply(cl = cl, X = lsl, FUN = function(x, m) {
    findmotif <- gsub("(^.*_)([[:alnum:]]{2})(\\.uniq.tab$)", "\\2", x, perl = TRUE)
    motif.in <- m[m$locus == findmotif, "motif"]
    system(sprintf("python microsatTabToseries.py -f %s -m %s", x, motif.in))
  }, m = motif)
}

if (do.chunk3) {
  ###############################################
  # CREATE CANDIDATE GENOTYPE FILE 
  ####
  # This part reads summarized "serie.tab" files from microsatTabToseries.py and produces candidate genotypes. 
  # Certain parts of it are run in parallel using the \code{parallel} package which comes shipped with R and 
  # works on Unix and Windows machines.
  # 
  # Input files should look like this (from microsatTabToseries.py):
  # id count sample.M1P5P.MM_065_P1
  # 1     K00209:87:HJKGNBBXX:6:1101:29153:1402_CONS_SUB_SUB  2727                    366
  # 2 K00209:87:HJKGNBBXX:6:1101:31061:1437_CONS_SUB_SUB_CMP  2423                    363
  # 3  K00209:87:HJKGNBBXX:6:1101:1986:1420_CONS_SUB_SUB_CMP   351                     84
  # 4     K00209:87:HJKGNBBXX:6:1101:23622:1455_CONS_SUB_SUB   292                     38
  # 5 K00209:87:HJKGNBBXX:6:1101:23886:1420_CONS_SUB_SUB_CMP    44                      1
  # 6      K00209:87:HJKGNBBXX:6:1101:4675:8541_CONS_SUB_SUB    44                      3
  # sample.M1P5P.MM_065_P2 sample.M1P5P.MM_065_P3 sample.M1P5P.MM_065_P4 sample.M1P5P.MM_065_P5
  # 1                    319                    269                    203                    449
  # 2                    270                    125                     85                    409
  # 3                     37                     16                     11                     47
  # 4                     50                     39                     16                     46
  # 5                      2                      1                      0                      4
  # 6                      5                      0                      0                      1
  # sample.M1P5P.MM_065_P6 sample.M1P5P.MM_065_P7 sample.M1P5P.MM_065_P8 seq_length
  # 1                    290                    291                    540         69
  # 2                    424                    482                    265         77
  # 3                     45                     68                     43         73
  # 4                     24                     24                     55         65
  # 5                      1                      4                     31         81
  # 6                     33                      1                      1         69
  # sequence series
  # 1             ttcttttcttttctttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      1
  # 2     ttcttttcttttctttctttctttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      1
  # 3         ttcttttcttttctttctttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      1
  # 4                 ttcttttcttttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      1
  # 5 ttcttttcttttctttctttctttctttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      1
  # 6             tccttttcttttctttctttctttctttctttctttctttctttctttctttctttctttcaacaaaca      2
  # 
  #  And output is:
  #       sample run count_run locus allele
  # 3  M0P0K.MM  P1         9   064   97_1
  # 9  M0P0K.MM  P2        88   064   97_1
  # 15 M0P0K.MM  P3         5   064   97_1
  # 21 M0P0K.MM  P4        51   064   97_1
  # 27 M0P0K.MM  P5        17   064   97_1
  # 33 M0P0K.MM  P6        92   064   97_1
  # 39 M0P0K.MM  P7        16   064   97_1
  # 45 M0P0K.MM  P8        66   064   97_1
  # sequence
  # 3  gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 9  gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 15 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 21 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 27 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 33 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 39 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  # 45 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
  ###############################################
  
  library(data.table)
  
  # fetch all series files and append some columns to aid in debugging and viewing of the data
  xy <- data.frame(files = list.files(dir.lsl, pattern = "serie.tab$", full.names = TRUE),
                   stringsAsFactors = FALSE)
  xy$lib <- gsub("^MICROSAT\\.PCR_(DAB\\d+)_.*$", "\\1", basename(xy$files))
  xy$sample <- gsub("^.*_(.*)_\\d{2}_serie\\.tab$", "\\1", xy$files)
  xy$locus <- gsub("^.*_(.*)_(\\d{2})_serie\\.tab$", "\\2", xy$files)
  
  smp <- unique(xy$sample)
  xy[xy$sample %in% sample(smp, size = 1), ]
  
  # weird samples
  ws <- sapply(split(xy, f = xy$sample), nrow)
  
  # Notice that some samples have been used in more than 1 library.
  message("Number of loci per sample.")
  print(table(ws))
  
  message("Samples with 156 files.")
  print(ws[ws == 156])
  
  message("Samples with 26 files.")
  print(ws[ws == 26])
  
  message("Samples with 39 files.")
  print(ws[ws == 39])
  
  message(sprintf("Processing %d files.", nrow(xy)))
  
  # For each file, read in the data, sort it according to total count of reads, remove unnecessary columns,
  # reflow ata into a long format for all repeats, modify or create some new columns like sample, locus, run,
  # and finish by reordering the result according to count within each run.
  
  xy$proc <- sprintf("%s (%d/%d)", xy$files, 1:nrow(xy), nrow(xy))
  xy <- split(xy, f = 1:nrow(xy))
  
  genotypes <- sapply(X = xy, FUN = function(x) {
    # Prints progress report every 100 files.
    if (as.numeric(gsub("^.*\\((\\d+)/\\d+\\)$", "\\1", x$proc)) %% 100 == 0) {
      message(sprintf("Processing %s at locus %s at %s", x$proc, x$locus, Sys.time()))
    }
    
    # Import data.
    io <- fread(x$files)
    
    # If there are no reads in the file, return NA.
    if (sum(io$count) == 0) {
      message(sprintf("No reads in %s.", x$files))
      return(NA)
    }
    
    # Subset six most common alleles (this can be varied).
    io <- io[order(io$count, decreasing = TRUE), ]
    
    # Remove unnecessary columns.
    io$id <- NULL # uuid to be added later
    io$count <- NULL # sum of count_run ~ seq_length + series
    
    # Reflow data into a long format.
    io <- gather(io, key = run, value = count_run, 
                 c(-sequence, -seq_length, -series))
    
    # Add sample names, locus, allele and run number.
    # Explanation of regex:
    # string should begin with "sample.", then find a *group*, then find _,
    # then find digits, followed by _, then by P, followed by more digits
    # and finally the string should end on that. Extract *group*.
    # The next two regex statements are similar, except they vary in number
    # of position of groups. \\d{1} means "find 1 digit".
    io$lib <- gsub("^MICROSAT\\.PCR_(.*)_(.*)_(\\d+)_serie.tab$", "\\1", basename(x$files))
    io$sample <- gsub("^MICROSAT\\.PCR_(.*)_(.*)_(\\d+)_serie.tab$", "\\2", basename(x$files))
    io$locus <- gsub("^MICROSAT\\.PCR_(.*)_(.*)_(\\d+)_serie.tab$", "\\3", basename(x$files))
    io$position <- gsub("^sample:(.*)_(\\d+)_(PP\\d{1})", "\\2", io$run)
    io$run <- gsub("^sample:(.*)_(\\d+)_(PP\\d{1})", "\\3", io$run)
    io$allele <- paste(io$seq_length, io$series, sep = "_")
    
    # Notice that `serie` and `seq_length` are no longer included as they are
    # encapsulated in `allele`.
    io <- io[, c("sample", "run", "count_run", "locus", "allele", "sequence", "lib", "position")]
    
    # Sort by run and then decreasing count_run.
    io <- io[order(rev(io$run), io$count_run, decreasing = TRUE), ]
    io <- io[io$count_run > 0, ]
    io
  }, simplify = FALSE)
  
  ## Save data into a .RData file in case shit hits the fan.
  genotypes <- genotypes[!is.na(genotypes)]
  save(genotypes, file = "./DAB/data/raw_genotypes_dab_hiseq1.RData")
}

if (do.chunk4) {
  #### Prepare data for saving. ####
  ##################################
  require(data.table)
  load("./DAB/data/raw_genotypes_dab_hiseq1.RData")
  genotypes <- rbindlist(genotypes)
  
  # Save blk data into its own data.frame for tidy purposes.
  # genotypes$sequence <- as.character(genotypes$sequence)
  gt <- genotypes # I feel typing "gt" is faster than "genotypes"
  rm(genotypes)
  
  # check if each sequence has only one locus
  # unique(aggregate(locus ~ sequence, data = gt, FUN = function(x) length(unique(x)))$locus)
  # ... expecting [1] 1
  
  # Standardize loci names
  gt.by <- by(data = gt, INDICES = list(gt$locus), FUN = function(x) {
    x <- x[order(x$count_run, decreasing = TRUE), ]
    x$new_allele <- as.numeric(factor(x$sequence, levels = unique(x$sequence), 
                                      labels = 1:length(unique(x$sequence))))
    x$new_allele <- paste(gsub("^(\\d+)_\\d$", "\\1", x$allele), x$new_allele, sep = "_")
    x[order(x$sample, x$run, rev(x$count_run)), ]
  })
  
  gt <- rbindlist(gt.by)
  # rm(gt.by)
  setDT(gt)
  
  # feel free to filter out junk
  gt <- gt[gt$count_run > 3, ]
  
  gt[, length := as.numeric(gsub("^(\\d+)_.*$", "\\1", gt$new_allele))]
  
  getcols <- c("sample", "run", "count_run", "locus", "lib", "length", "position", "sequence")
  gt <-  gt[, ..getcols]
  
  gt <- gt[order(gt$sample, gt$locus, gt$run, rev(gt$count_run))]
  gt$run <- gsub("PP", "", gt$run) # remove PP for prettier printing
  
  # rename columns
  names(gt) <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", "length", "Position", "Sequence")
  
  # add tag combo
  xy <- sapply(list.files("./DAB/1_ngsfilters_hiseq1/", pattern = ".ngsfilter", full.names = TRUE),
               FUN = read.table, simplify = FALSE)
  xy <- do.call(rbind, xy)
  rownames(xy) <- NULL
  xy <- as.data.table(xy[, c(1, 2, 3)])
  names(xy) <- c("V1", "V2", "TagCombo")
  
  
  gt[, fn := sprintf("%s_%s_PP%s", Sample_Name, Position, Plate)]
  gt[, fl := sprintf(sprintf("UA_MxRout1_%s", Marker))]
  
  gt <- merge(gt, xy, by.x = c("fn", "fl"), by.y = c("V2", "V1"))
  gt[, fn := NULL]
  gt[, fl := NULL]
  
  save(gt, file = "./DAB/data/genotypes_dab_hiseq1_cleaned.RData")
}

if (do.chunk5) {
  # Clean genotypes ####
  #####
  require(data.table)
  library(fishbone)
  
  load("./DAB/data/genotypes_dab_hiseq1_cleaned.RData")
  
  mt <- fread("./DAB/pars.csv", dec = ",",
              colClasses = list(character = c(1, 2), numeric = c(3, 4, 5, 6)),
              stringsAsFactors = FALSE, header = TRUE)
  gen <- split(gt, f = list(gt$Sample_Name, gt$Marker, gt$Plate))
  
  empty <- sapply(gen, nrow)
  gen.empty <- names(gen[empty])
  gen[empty == 0] <- NULL
  gen2 <- sapply(gen[1:100], FUN = callAllele, tbase = mt, simplify = FALSE)
}
