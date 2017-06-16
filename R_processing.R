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
  inputfile <- list.files(pattern = "MICROSAT")
  
  # extract sample names
  samplenames <- as.character(unique(sn$V2))
  samplenames <- unique(sapply(strsplit(samplenames, "_"), "[", 1))
  
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
    microsatTabExtract(filename = x["inputfile"] , sampleName = x["samplename"])
  })
  
  # apply(X = sample.loci.combo, MARGIN = 1, FUN = function(x) {
  #   microsatTabExtract(filename = x["inputfile"] , sampleName = x["samplename"])
  # })
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
  #########################################
  #### Step 3: find consensus genotype ####
  #########################################
  
  # This script select for each repeat up to six potential alleles and print the data into a table.
  # It requires 5 arguments : 1) the number of repeats for each sample 2) the microsat motif length for this marker 
  # 3) the number of samples 4) the minimal number of sequences per serie for a sample 5) the name of the input tab 
  # file (output from microsatTabToserie.py)	
  # 
  # Ex : 
  #   R -q --vanilla --args 8 4 19 10 MICROSAT.PCR_UARef_03_serie.tab < microsat_alleles.r
  #                         | | |  |__ minimal number of sequences per series for a sample
  #                         | | |_____ number of samples
  #                         | |_______ microsat motif length for this marker (e.g. 03)
  #                         |_________ number of repeats for each sample
  #
  #
  # It will generate the file : table2n_MICROSAT.PCR_UARef_03_serie.tab
  
  message("Performing step 3.")
  
  oldwd <- getwd()
  setwd("./data")
  
  xys <- list.files(pattern = "^MICROSAT\\.PCR.*_serie.tab")
  message(sprintf("Found %s files to be processed", length(xys)))
  print(head(xys))
  
  # its <- grepl("PCR_POS_14", xys)
  # sapply(X = xys[500], FUN = function(x, m, ns, nr) {
  #   message(sprintf("Processing file %s", x))
  #   tm <- gsub("(^MICROSAT\\.PCR)_([[:alnum:]\\.]*)_([[:alnum:]]{2})(\\_serie\\.tab$)", "\\3", x, perl = TRUE)
  #   ml <- nchar(m[m$locus == tm, "motif"])
  #   
  #   microsatAlleles(nrep = nr, motifLength = ml, numsamples = ns,
  #                   numseqfamily = 10, input = x)
  # }, m = motif, ns = num.samples, nr = num.repeats)
  
  parSapply(cl = cl, X = xys, FUN = function(x, m, ns, nr) {
    tm <- gsub("(^MICROSAT\\.PCR)_([[:alnum:]\\.]*)_([[:alnum:]]{2})(\\_serie\\.tab$)", "\\3", x, perl = TRUE)
    ml <- nchar(m[m$locus == tm, "motif"])
    
    microsatAlleles(nrep = nr, motifLength = ml, numsamples = ns,
                    numseqfamily = 10, input = x)
  }, m = motif, ns = num.samples, nr = num.repeats)
}
if (FALSE) {
  #################
  #### Step 4: ####
  #################
  
  # This script defines a consensus genotype, if possible, using data from the different repeats for a sample.
  # It requires 5 arguments : 1) the number of repeats for each sample 2) the number of samples 3) the threshold 
  # to apply to seq counts 4) the homozygous threshold 5) the name of the input file (output from microsat_alleles.r)
  # The homozygous threshold corresponds to the percentage of repeats matching an homozygous genotype that are
  # required to validate an homozygous genotype). In the following example, this threshold is setup at 0.5, so at 
  # least 4 out of the 8 repeats are required to validate an homozygous genotype.
  # 
  # 
  # Ex:
  #   R –q --vanilla --args 8 19 10 0.5 table2n_MICROSAT.PCR_UARef_03_serie.tab < microsat_consensus.r
  #                         | |  |   |_ homozygous threshold
  #                         | |  |_____ threshold for sequence count
  #                         | |________ number of samples
  #                         |__________ number of repeats
  #                               
  # It will generate the files : ctable2n_MICROSAT.PCR_UARef_03_serie.tab, 
  # cns_table2n_MICROSAT.PCR_UARef_03_serie.tab, cnsh_table2n_MICROSAT.PCR_UARef_03_serie.tab and UA_3_ref 
  # (if this latter was not provided)
  # 
  # A reference allele file can be provided for allele names to keep identical names for the alleles from one 
  # experiment to the other. The program first output the consensus genotypes from the analyzed run 
  # (ctable2n_MICROSAT.PCR_UARef_03_serie.tab  and cns_table2n_MICROSAT.PCR_UARef_03_serie.tab) and then take 
  # in account the reference allele file of the marker to change allele names according to this file if required
  # and output the results in a new file (cnsh_table2n_MICROSAT.PCR_UARef_03_serie.tab).
  # If no reference allele file is provided, the program creates a new one and print into this file the alleles 
  # found in this run (in this case cns_ and cnsh_ files are equal).
  # If a reference allele file is provided and new alleles found in the run, the reference allele file is updated.
  
  # oldwd <- getwd()
  # setwd("./data")
  
  xyt <- list.files(pattern = "^table2n_")
    sapply(xyt, FUN = microsatConsensus, repeats = num.repeats, 
           NbEch = num.samples,
           threshold = 10, 
           homozygousThreshold = 0.5)
}

if (FALSE) {
  # Creating the final output.
  xy <- list.files(pattern = "ctable2n")
  
  gts <- sapply(xy, FUN = extractGenotypes, simplify = FALSE)
  gts.orig <- gts
  gts <- do.call(rbind, gts)
  
  write.table(gts, file = "genotypes_UA_GATC.txt", sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE)
}
