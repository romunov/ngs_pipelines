library(parallel)
source("microsatAlleles.R")
source("microsatConsensus.R")
source("microsatTabExtract.R")

# Initial chunk, if () used to fold unnecessary code.
if (TRUE) {
# start cores
  if (Sys.info()["sysname"] == "Windows") {
    ncores <- 4
  } else {
    ncores <- 46
  }
  
  message(sprintf("Powering up %d cores.", ncores))
  cl <- makeCluster(ncores)
  
  clusterExport(cl = cl, varlist = c("microsatAlleles", "microsatTabExtract"))
  on.exit(stopCluster(cl))
  
  # load data and sample names
  message("Importing files.")
  sn <- list.files(pattern = ".ngsfilter", full.names = TRUE)
  sn <- read.table(sn)
  inputfile <- list.files("./data", pattern = "^MICROSAT\\.PCR_UA_.*\\.tab$", full.names = TRUE)
  
  # extract sample names
  samplenames <- as.character(unique(sn$V2))
  samplenames <- gsub("^(.*)_[[:alnum:]]{3}_P\\d+$", "\\1", samplenames)
  
  # load motif data
  motif <- read.table("locus_motifs.txt", header = TRUE,
                      colClasses = c("character", "character"))
  message("Files imported.")
  
  # Script levels parameters:
  num.samples <- length(unique(samplenames))
  num.repeats <- 8
}

if (FALSE) {
  ############################################
  #### Step 1: Create sample-locus files. ####
  ############################################
  
  # Example:
  # R -q --vanilla --args MICROSAT.PCR_UA_MplexRout1_03.tab UARef < microsat_tabextract.r
  #                         |_ locus file with all samples    |_ sample name
  
  # to extract data for all loci of all samples, we need to create all
  # combinations
  message("Performing step 1.")
  
  sample.loci.combo <- expand.grid(samplename = samplenames, inputfile = inputfile)
  
  print(head(sample.loci.combo))
  
  parApply(cl = cl, X = sample.loci.combo, MARGIN = 1, FUN = function(x) {
    microsatTabExtract(filename = x["inputfile"] , samplename = x["samplename"])
  })
  
  apply(X = sample.loci.combo[5000, , drop = FALSE], MARGIN = 1, FUN = function(x) {
    microsatTabExtract(filename = x["inputfile"] , samplename = x["samplename"])
  })
}
if (FALSE) {
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
  
  # find files created by microsat_tabextract.r
  # NOTICE: at this point, ZF locus is excluded from the analysis
  samplelocus <- list.files(pattern = "MICROSAT.PCR_[^UA].*_\\d{2}\\.uniq.tab$")
  
  message(sprintf("Found %d files to process.", length(samplelocus)))
  
  parSapply(cl = cl, X = samplelocus, FUN = function(x, m) {
    findmotif <- gsub("(^.*_)([[:alnum:]]{2})(\\.uniq.tab$)", "\\2", x, perl = TRUE)
    motif.in <- m[m$locus == findmotif, "motif"]
    
    # message(sprintf("Processing %s", x))
    # shell(sprintf("microsatTabToseries.py -f %s -m %s", x, motif.in))
    system(sprintf("python microsatTabToseries.py -f %s -m %s", x, motif.in))
  }, m = motif)
  
  # sapply(X = samplelocus, FUN = function(x, m) {
  #   findmotif <- gsub("(^.*_)([[:alnum:]]{2})(\\.uniq.tab$)", "\\2", x, perl = TRUE)
  #   motif.in <- m[m$locus == findmotif, "motif"]
  #   
  #   message(sprintf("Processing %s", x))
  #   # shell(sprintf("microsatTabToseries.py -f %s -m %s", x, motif.in))
  #   system(sprintf("python microsatTabToseries.py -f %s -m %s", x, motif.in))
  # }, m = motif)
}

if (FALSE) {
  ##################################
  # CREATE CANDIDATE GENOTYPE FILE #
  ##################################
  # Make sure to adapt the creation of `xy` to get the correct files.
  source("analysis_get_genotypes.R")
}
