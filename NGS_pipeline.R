# This is the R part of the pipeline used to analyze NGS data. Study
# species is brown bear, Ursus arctos.

library(parallel)
library(tidyr)
library(data.table)
source("custom_functions.R")
source("determine_sex_from_ZF.R")  # needed to call sex from samples

# Have these files in the "home" folder of the analysis (ngs_pipelines repo)
source("microsatTabExtract.R")
locus.motifs <- "locus_motifs.txt"  # motif represented by each locus
microsattotab.py <- "microsatTabToseries.py"  # a python script which finds "series"

do.chunk.init <- TRUE  # this will be set to TRUE if needed by specific chunk
do.chunk1 <- TRUE  # create files for each sample/locus using microsatTabExtract
do.chunk2 <- TRUE  # find series using python script
do.chunk3 <- TRUE  # prepare candidate alleles
do.chunk4 <- TRUE  # clean candidate alleles
do.chunk5 <- TRUE  # call alleles
do.chunk6 <- TRUE  # call sex

# for chunks running in parallel, always turn on init chunk to start up workers
if (do.chunk1) do.chunk.init <- TRUE
if (do.chunk2) do.chunk.init <- TRUE

#### USER INPUT ####
# Specify project name.
proj.name <- "FCC2"
# Specify output file which will be placed inside /data of the project folder.
raw.rdata <- sprintf("./%s/data/raw_genotypes_%s.RData", proj.name, proj.name)
raw.cleaned.rdata <- sprintf("./%s/data/genotypes_%s_cleaned.RData", proj.name, proj.name)
raw.final <- sprintf("./%s/data/final_dab_%s.RData", proj.name, proj.name)
raw.final.txt <- sprintf("./%s/data/%s_genotypes.txt", proj.name, proj.name)

# specify (intermediate() folder names.
dir.ngsfilter <- sprintf("./%s/1_ngsfilters", proj.name)
dir.uniq.tab <- sprintf("./%s/2_uniq_tab", proj.name)
dir.lsl       <- sprintf("./%s/3_lib_sample_locus", proj.name)

#### PROCESSING OF USER INPUT ####
# If project folder doesn't exist yet, create on.
if (!dir.exists(proj.name)) {
  dir.create(proj.name)
  dir.create(sprintf("./%s/data", proj.name))
} else {
  # create /data inside the project name
  if (!dir.exists(sprintf("./%s/data", proj.name))) {
    dir.create(sprintf("./%s/data", proj.name))
  }
}


#### INITIALIZE PARALLEL TOOLS ####
if (do.chunk.init) {
  ncores <- detectCores() - 1
  
  message(sprintf("Powering up %d processes.", ncores))
  cl <- makeCluster(ncores, outfile = "clusterfuck.txt")
  on.exit(stopCluster(cl))
}

if (!dir.exists(dir.ngsfilter)) {
  dir.create(dir.ngsfilter)
}

if (!dir.exists(dir.uniq.tab)) {
  dir.create(dir.uniq.tab)
}

if (!dir.exists(dir.lsl)) { # create if doesn't exist
  dir.create(dir.lsl)
}

#### BEGIN PROCESSING CHUNKS ####
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
  sn <- sapply(sn, fread, header = FALSE, colClasses = "character", simplify = FALSE)
  
  # TODO: this part could be programmed to move .uniq.tab files from .. to ./PROJECT/2_uniq_tab.
  inputfile <- list.files(dir.uniq.tab, pattern = "^MICROSAT.*\\.uniq\\.tab$", full.names = TRUE)
  libnum <- gsub("^.*-(\\d+)_UA_.*\\.uniq.tab$", "\\1", basename(inputfile))
  # libnum <- sprintf("%02d", as.numeric(libnum))  # this line is needed if lib number is not preceded by 0
  
  inputfile <- split(inputfile, f = libnum)
  
  # Remove any ngs filters for which there is no data (e.g.).
  for (i in 1:length(sn)) {
    # Check if any of the library numbers appear in the ngsfilter name.
    lib.detected <- any(sapply(unique(libnum), grepl, x = names(sn)[i]))
    # If library is not detected, remove it from ngs filters
    if (!lib.detected) {
      message(sprintf("Library %s not detected, excluded from the analysis.", names(sn)[i]))
      sn[i] <- NULL
    }
  }
  # check that correct ngs filter (sample names) are assigned to correct library data
  message("Do input numbers match library ngsfilter file?")
  print(data.table(ngsfilter = names(sn), input = names(inputfile)))
  
  # extract sample names
  # TODO: sample names should be tied to tag combo or position if we want to prevent pooling
  # of repeats between plates of a given run
  samplenames <- sapply(sn, FUN = function(x) as.character(x$V2), simplify = FALSE)
  # samplenames <- sapply(samplenames, FUN = function(x) {
  #   unique(gsub("^(.*)_[[:alnum:]]{3}_PP\\d+$", "\\1", x))
  # }, simplify = FALSE)
  
  # Create a combination to extract alleles from all loci for a given sample.
  slc <- mapply(FUN = function(x.sample, x.input) {
    expand.grid(samplename = x.sample, inputfile = x.input)
  }, samplenames, inputfile, SIMPLIFY = FALSE)
  slc <- rbindlist(slc)
  slc$samplename <- as.character(slc$samplename)
  slc$inputfile <- as.character(slc$inputfile)
  
  message(sprintf("Processing %d files to %s.", nrow(slc), dir.lsl))
  
  # # Split the data into ncores chunks which are to be processed in parallel.  
  slc <- split(slc, f = slc$samplename)
  
  # progress report
  numsam <- length(slc)
  i <- 0
  
  out <- sapply(X = slc, FUN = function(x, dir.lsl, outdir) {
    i <<- i + 1
    message(sprintf("Processing %s (%d/%d)", unique(x[, samplename]), i, numsam)) # scoping out of current env!
    out <- parApply(cl = cl, X = x, MARGIN = 1, FUN = function(y, outdir) {
      # out <- apply(X = x, MARGIN = 1, FUN = function(y, outdir) {
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
  # lsl = library-sample-locus files
  
  oldwd <- getwd()
  # NOTICE: copy microsatTabToseries.py to where the uniq.tab files are being processed
  file.copy(from = microsattotab.py, to = paste(dir.lsl, "microsatTabToseries.py", sep = "/"),
            overwrite = TRUE)
  
  lsl <- list.files(path = dir.lsl, pattern = "_\\d+\\.uniq\\.tab$")
  # load motif data
  motif <- read.table(locus.motifs, header = TRUE,
                      colClasses = c("character", "character"))
  
  message(sprintf("Found %d files to process.", length(lsl)))
  message("Sample of files:")
  print(head(lsl))
  
  if (!file.exists("microsatTabToseries.py")) stop("microsatTabToseries.py not found")
  
  # nodes for some reason are not set with current working dir, so do it explicitly
  # on each node
  clusterExport(cl = cl, varlist = "dir.lsl")
  clusterEvalQ(cl = cl, expr = {setwd(dir.lsl)})
  
  setwd(dir.lsl)
  parSapply(cl = cl, X = lsl, FUN = function(x, m) {
    # sapply(X = lsl, FUN = function(x, m) {
    findmotif <- gsub("(^.*_)([[:alnum:]]{2})(\\.uniq.tab$)", "\\2", x, perl = TRUE)
    motif.in <- m[m$locus == findmotif, "motif"]
    system(sprintf("python microsatTabToseries.py -f %s -m %s", x, motif.in))
  }, m = motif)
  
  # This chunk is used for debugging.
  # sapply(X = lsl, FUN = function(x, m) {
  #   findmotif <- gsub("(^.*_)([[:alnum:]]{2})(\\.uniq.tab$)", "\\2", x, perl = TRUE)
  #   motif.in <- m[m$locus == findmotif, "motif"]
  #   system(sprintf("python microsatTabToseries.py -f %s -m %s", x, motif.in))
  # }, m = motif)
  
  setwd(oldwd)
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
  #        sample                       run count_run locus allele
  # 1  B07blank01 sample:B07blank01_072_PP1         9    03   59_1
  # 2  B07blank01 sample:B07blank01_072_PP1         8    03   63_1
  # 3  B07blank01 sample:B07blank01_072_PP1         3    03   63_2
  # 4  B07blank01 sample:B07blank01_072_PP1         2    03   67_1
  # 5  B07blank01 sample:B07blank01_072_PP1         2    03   59_3
  # 6  B07blank01 sample:B07blank01_072_PP1         1    03   67_4
  # 7  B07blank01 sample:B07blank01_072_PP1         1    03   63_5
  # 8  B07blank01 sample:B07blank01_072_PP1         1    03   63_6
  # 9  B07blank01 sample:B07blank01_072_PP1         1    03   59_7
  # 10 B07blank01 sample:B07blank01_072_PP1         1    03   59_8
  # 11 B07blank01 sample:B07blank01_072_PP1         1    03   55_1
  #                                                               sequence    lib
  # 1          aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctc DAB045
  # 2      aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctc DAB045
  # 3      aaatcctgtaacaaatctatctatctatctatctatctatctatctacctatctatctatctc DAB045
  # 4  aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctatctc DAB045
  # 5          aaatcctgtgacaaatctatctatctatctatctatctatctatctatctatctatctc DAB045
  # 6  aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctatcta DAB045
  # 7      aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatccatctatctc DAB045
  # 8      aaatcctgtaacaaatctatctatctatctatctatctatccatctatctatctatctatcta DAB045
  # 9          aaatcctgtaacgaatctatctatctatctatctatctatctatctatctatctatctc DAB045
  # 10         aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatcac DAB045
  # 11             aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctc DAB045
  #    position plate
  # 1       072   PP1
  # 2       072   PP1
  # 3       072   PP1
  # 4       072   PP1
  # 5       072   PP1
  # 6       072   PP1
  # 7       072   PP1
  # 8       072   PP1
  # 9       072   PP1
  # 10      072   PP1
  # 11      072   PP1
  ###############################################
  
  # fetch all series files and append some columns to aid in debugging and viewing of the data
  # (!) notice that this step does not account for locus ZF
  stringy <- "^MICROSAT.PCR_(DAB\\d+)_(.*)_(.*)_(.*)_(.*)_serie.tab$"
  #                                                    |- locus    (5)
  #                                               |------ plate    (4)
  #                                          |----------- position (3)
  #                                     |---------------- sample   (2)
  #                             |------------------------ library  (1)
  
  xy <- data.table(files = list.files(dir.lsl, pattern = "serie.tab$", full.names = TRUE),
                   stringsAsFactors = FALSE)
  xy$lib <- gsub(stringy, "\\1", basename(xy$files))
  xy$sample <- gsub(stringy, "\\2", basename(xy$files))
  
  smp <- unique(xy$sample)
  xy[xy$sample %in% sample(smp, size = 1), ]
  
  # weird samples
  ws <- sapply(split(xy, f = xy$sample), nrow)
  tws <- table(ws)
  
  # Notice that some samples have been used in more than 1 library.
  message("Number of loci per sample.")
  sapply(names(tws)[2:length(tws)], FUN = function(x) {
    ws[ws == x]
  })
  
  message(sprintf("Processing %d files.", nrow(xy)))
  
  # For each file, read in the data, sort it according to total count of reads, remove unnecessary columns,
  # reflow data into a long format for all repeats, modify or create some new columns like sample, locus, run,
  # and finish by reordering the result according to count within each run.
  
  xy$proc <- sprintf("%s (%d/%d)", xy$files, 1:nrow(xy), nrow(xy))
  xy <- split(xy, f = 1:nrow(xy))
  
  genotypes <- sapply(X = xy, FUN = function(x, ptn) {
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
    
    # Add library, sample names, locus, position, allele and run number.
    # Explanation of regex:
    # string should begin with "sample.", then find a *group*, then find _,
    # then find digits, followed by _, then by P, followed by more digits
    # and finally the string should end on that. Extract *group*.
    # The next two regex statements are similar, except they vary in number
    # of position of groups. \\d{1} means "find 1 digit".
    
    io$lib <- gsub(stringy, "\\1", basename(x$files))
    io$sample <- gsub(stringy, "\\2", basename(x$files))
    io$locus <- gsub(stringy, "\\5", basename(x$files))
    io$position <- gsub(stringy, "\\3", basename(x$files))
    io$run <- gsub(stringy, "\\4", basename(x$files))
    io$locus <- gsub(stringy, "\\5", basename(x$files))
    io$allele <- paste(io$seq_length, io$series, sep = "_")
    
    # Notice that `serie` and `seq_length` are no longer included as they are
    # encapsulated in `allele`.
    io <- io[, c("sample", "run", "count_run", "locus", "allele", "sequence", "lib", "position")]
    
    # Sort by run and then decreasing count_run.
    io <- io[order(rev(io$run), io$count_run, decreasing = TRUE), ]
    io <- io[io$count_run > 0, ]
    io
  }, ptn = stringy, simplify = FALSE)
  
  ## Save data into a .RData file in case shit hits the fan.
  genotypes <- genotypes[!is.na(genotypes)]
  save(genotypes, file = raw.rdata)
}

if (do.chunk4) {
  
  message(sprintf("(%s) Processing chunk #4", Sys.time()))
  ##################################
  #### Prepare data for saving. ####
  ##################################
  if (!do.chunk3) {
    require(data.table)
    load(raw.rdata)
  }
  gt <- rbindlist(genotypes)
  
  # # uncomment only if lib, position and run are FUBAR
  # gt$lib <- gsub("(^DAB\\d+)_A\\d+$", "\\1", gt$lib)
  # gt$position <- gsub("^sample:.*_(\\d{3})_PP\\d$", "\\1", gt$position)
  # gt$run <- gsub("^sample:.*_(PP\\d)$", "\\1", gt$run)
  
  # Standardize loci names
  gt.by <- by(data = gt, INDICES = list(gt$locus), FUN = function(x) {
    x <- x[order(x$count_run, decreasing = TRUE), ]
    x$new_allele <- as.numeric(factor(x$sequence, levels = unique(x$sequence),
                                      labels = 1:length(unique(x$sequence))))
    x$new_allele <- paste(gsub("^(\\d+)_\\d$", "\\1", x$allele), x$new_allele, sep = "_")
    x[order(x$sample, x$run, rev(x$count_run)), ]
  })
  gt <- rbindlist(gt.by)
  rm(gt.by)
  
  # feel free to filter out junk
  gt <- gt[gt$count_run > 9, ]
  
  gt[, length := as.numeric(gsub("^(\\d+)_.*$", "\\1", gt$new_allele))]
  
  # subset only relevant columns
  getcols <- c("sample", "run", "count_run", "locus", "lib", "length", "position", "sequence")
  gt <-  gt[, ..getcols]
  
  gt <- gt[order(gt$sample, gt$locus, gt$run, rev(gt$count_run))]
  gt$run <- gsub("^.*(\\d+)$", "\\1", gt$run) # extract run number for prettier printing
  
  # rename columns
  names(gt) <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", "length", "Position", "Sequence")
  
  # add tag combo
  xy <- sapply(list.files(dir.ngsfilter, pattern = ".ngsfilter", full.names = TRUE),
               FUN = function(x) {
                 out <- fread(x, header = FALSE, colClasses = "character")
                 out$rn <- gsub("\\.ngsfilter", "", basename(x))
                 out
               },
               simplify = FALSE)
  
  # Remove any ngs filters for which there is no data (e.g.).
  for (i in 1:length(xy)) {
    # Check if any of the library numbers appear in the ngsfilter name.
    lib.detected <- any(sapply(unique(libnum), grepl, x = names(xy)[i]))
    # If library is not detected, remove it from ngs filters
    if (!lib.detected) {
      message(sprintf("Library %s not detected, excluded from the analysis.", names(xy)[i]))
      xy[i] <- NULL
    }
  }
  
  xy <- rbindlist(xy)
  rownames(xy) <- NULL
  xy <- as.data.table(xy[, c(1, 2, 3, 8)])
  names(xy) <- c("V1", "V2", "TagCombo", "rn")
  
  # create columns by which to merge
  gt[, fn := sprintf("%s_%s_PP%s", Sample_Name, Position, Plate)]
  # this is here temorarily - switch between primer plate designation P and PP
  # gt[, fn := sprintf("%s_%s_P%s", Sample_Name, Position, Plate)]
  gt[, fl := sprintf(sprintf("UA_MxRout1_%s", Marker))]
  
  gt <- merge(gt, xy, by.x = c("fn", "fl", "Run_Name"), by.y = c("V2", "V1", "rn"))
  gt[, fn := NULL]
  gt[, fl := NULL]
  
  message(sprintf("(%s) Chunk4: Writing data to file.", Sys.time()))
  save(gt, file = raw.cleaned.rdata)
  message(sprintf("(%s) Done processing chunk #4", Sys.time()))
}

if (do.chunk5) {
  #### Clean genotypes and calling alleles ####
  # Accepts data in the following format:
  
  # Sample_Name Plate Read_Count Marker Run_Name length Position
  # 1:       Blk01     1          9     03    DAB01     55      040
  # 2:       Blk01     1          5     03    DAB01     59      040
  # 3:       Blk01     1          4     03    DAB01     67      040
  # 4:       Blk01     1         10     16    DAB01     79      040
  # 5:       Blk01     1          9     16    DAB01     75      040
  # ---                                                             
  #   23361:    M20T8.MM     8          6     63    DAB01     84      069
  # 23362:    M20T8.MM     8          5     63    DAB01     77      069
  # 23363:    M20T8.MM     8         37     65    DAB01     93      069
  # 23364:    M20T8.MM     8         18     65    DAB01     89      069
  # 23365:    M20T8.MM     8         11     65    DAB01     81      069
  # Sequence          TagCombo
  # 1:                                       aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctc tagtcgca:gatcgcga
  # 2:                                   aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctc tagtcgca:gatcgcga
  # 3:                           aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctatctc tagtcgca:gatcgcga
  # 4:               aattttgttctttctttctttctttctttctttctttctttctttctttctttctttctttctttctctttcttttcag tagtcgca:gatcgcga
  # 5:                   aattttgttctttctttctttctttctttctttctttctttctttctttctttctttctttctctttcttttcag tagtcgca:gatcgcga
  # ---                                                                                                                
  #   23361:          tccatccatcatccatcatccatccatccatccatccatccatccatccatccatccggttactgctcatttaaaagcatggtc ggatagca:gtgatctc
  # 23362:                 tccatccatcatccatccatccatccatccatccatccatccatccatccggttactgctcatttaaaagcatggtc ggatagca:gtgatctc
  # 23363: gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc ggatagca:gtgatctc
  # 23364:     gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc ggatagca:gtgatctc
  # 23365:             gaagcaacagggtatagatatatagagatagatagatagatagatagatagataaagagatttattataaggaattggctc ggatagca:gtgatctc
  #####
  message(sprintf("(%s) Chunk 5: Begin processing chunk.", Sys.time()))
  
  library(fishbone)
  
  if (!exists("gt")) load(raw.cleaned.rdata)
  
  data(mt) # from fishbone package
  system.time(out <- gt[, callAllele(c(.BY, .SD), tbase = mt, verbose = TRUE),
                        by = .(Sample_Name, Marker, Plate, Run_Name)])
  
  # data.table adds variables used to "by" - here we remove them
  out <- out[, 5:ncol(out)]
  
  message("Chunk 5: Writing final genotype to .RData file.")
  save(out, file = raw.final)
  
  message("Chunk 5: Writing final genotypes to a text file.")
  fwrite(out, file = raw.final.txt, sep = "\t")
  message("Done processing chunk #5.")
}

if (do.chunk6) {
  sex <- determineSex(proj.name = proj.name)
}
