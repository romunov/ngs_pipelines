library(data.table)
library(tools)

# This script will run the pipeline on all found libraries. You need to specify which file
# belongs to which library (see below for an example). During the pipeline, some intermediate
# files are created and can be discarded (or not). Some switches are implemented which enable
# you to throw them out. Beware that the size of the project will grow rapidly, partially on
# the account of these files.

# TODO: implement _moving_ of intermediate files to (a) predestined location(s)

#### Cleanup ####
remove.bad.illuminapairedend <- TRUE # remove files which contain Bad alignments
remove.unidentified.illuminapairended <- TRUE # remove files with unidentified reads
remove.pairended.illuminapairended <- TRUE # remove pairended result from illuminapairended
remove.alignement.illuminapairedend  <- TRUE # remove aligned files from illuminapairended
remove.filtered.files <- TRUE # after combining filtered files to a single file, remove split chunks?
####

# Find the number of available cores.
N <- as.integer(system("nproc", intern = TRUE))

# Read in data and split it according to library-run variable. This will take care of
# cases where one library has more repeats.
# File example:
#                                        file   library run                         filter
# 1: NG-12425_E45175_lib187741_5356_6_1.fastq DAB01GATC   1 dab_gatc_named_ctrls.ngsfilter
# 2: NG-12425_E45175_lib187741_5356_6_2.fastq DAB01GATC   1 dab_gatc_named_ctrls.ngsfilter
setwd("/home/romunov/Documents/ngs_gatc_06_2017")

xy <- fread("input_gatc.txt", header = TRUE)
xy <- split(xy, f = list(xy$library, xy$run))

# For each unique combination/library, run the following pipeline.
xy.split <- sapply(xy, FUN = function(x, N) {
  browser()
  file.no.ext <- file_path_sans_ext(x$file)
  # x is a data.frame with file, library and run columns
  
  # RUN #: run.splitfiles
  # split files into N parts
  runthis <- sprintf("fastqutils split %s %s %d", 
                     x$file,
                     file.no.ext,
                     N)
  # TODO: system(runthis)
  
  # RUN #: run.illuminapairedend
  # collect files
  strand2 <- file.no.ext[grepl("_2$", file.no.ext)]
  strand1 <- file.no.ext[grepl("_1$", file.no.ext)]
  
  chunks.strand1 <- list.files(pattern = "\\.\\d+\\.fastq$")
  chunks.strand1 <- data.frame(strand1 = chunks.strand1[grepl(strand1, chunks.strand1)])
  chunks.strand1$chunk <- gsub("^.*_1\\.(\\d+)\\.fastq$", "\\1", chunks.strand1$strand1)
  
  chunks.strand2 <- list.files(pattern = "\\.\\d+\\.fastq$")
  chunks.strand2 <- data.frame(strand2 = chunks.strand2[grepl(strand2, chunks.strand2)])
  chunks.strand2$chunk <- gsub("^.*_2\\.(\\d+)\\.fastq$", "\\1", chunks.strand2$strand2)
  
  chunks <- merge(chunks.strand2, chunks.strand1, by = "chunk")
  # these names will be used to pipe the result of illuminapairedend to
  chunks$output <- sprintf("rawdata_%s_%s.", unique(x$library), chunks$chunk)
  chunks$pairended <- sprintf("pairended_%s_%s.fastq", unique(x$library), chunks$chunk)
  
  # construct a string that will pairend forward and reverse read
  run.pe <- sprintf("illuminapairended -r %s %s | tee %s | obiannotate -S goodAli:'\"Alignement\" if score>40.00 else \"Bad\"' | obisplit -t goodAli -p %s",
                    # illuminapairedend
                    chunks$strand2, # forward strand
                    chunks$strand1, # reverse strand
                    # tee
                    # unique(x$library), # library
                    chunks$pairended, # which chunk
                    # obisplit
                    # unique(x$library), # library
                    chunks$output)
  
  # TODO: system(paste(run.pe, collapse = " & "))
  
  if (remove.bad.illuminapairedend) {
    rm.bad <- list.files(pattern = "^.*\\.Bad.fastq$")
    if (length(rm.bad) > 0) {
      file.remove(rm.bad)
    }
  }
  
  # Collect chunks which have been created by illuminapairedend (IPE) tools.
  chunkstring <- "^(rawdata_.*)_(\\d+)\\.Alignement\\.fastq$"
  get.ipe.chunks <- data.frame(input = list.files(pattern = "rawdata_.*.Alignement.fastq$"))
  get.ipe.chunks$chunk <- gsub(chunkstring, "\\2", get.ipe.chunks$input)
  get.ipe.chunks$root <- gsub(chunkstring, "\\1", get.ipe.chunks$input)
  get.ipe.chunks$filtered <- with(get.ipe.chunks, sprintf("%s_%s_filtered.fastq", root, chunk))
  get.ipe.chunks$unidentified <- with(get.ipe.chunks, sprintf("%s_%s_unidentified.fastq", root, chunk))
  
  # TODO: system(paste(sprintf("touch %s", get.ipe.chunks$input), collapse = " & "))
  
  # run .ngsfilter for each chunk
  ngsfilter <- sprintf("ngsfilter -t %s -u %s %s > %s",
                       unique(x$filter), # ngsfilter
                       get.ipe.chunks$unidentified, # unidentified
                       get.ipe.chunks$input, # input
                       get.ipe.chunks$filtered # output chunk number and
  )
  
  # TODO: system(paste(ngsfilter, collapse = " & "))
  
  # Remove files that hold unidentified sequences.
  # TODO: user should be also able to move the data to a prespecified location(s)
  if (remove.unidentified.illuminapairended) {
    rm.unident <- list.files(pattern = "^.*_unidentified.fastq$")
    if (length(rm.unident) > 0) {
      file.remove(rm.unident)
    }
  }
  
  if (remove.pairended.illuminapairended) {
    rm.parend <- list.files(pattern = "^rawdata_pairended.*\\.fastq$")
    if (length(rm.parend) > 0) {
      file.remove(rm.parend)
    }
  }
  
  if (remove.alignement.illuminapairedend) {
    rm.alg <- list.files(pattern = "^rawdata_.*Alignement\\.fastq$")
    if (length(rm.alg) > 0) {
      file.remove(rm.alg)
    }
  }
  
  # merge chunks into a single chunk before processing further
  filtered <- list.files(pattern = "^rawdata.*_filtered\\.fastq$")
  # system(sprintf("cat %s > %s_filtered.fastq", 
  #                paste(filtered, collapse = " "),
  #                unique(x$library))
  # )
  
  if (remove.filtered.files) {
    if (length(filtered) > 0) {
      file.remove(filtered)
    }
  }
  
  # RUN #: run.splitdatabylocus
  # split data by locus
  run.splitlocus <- sprintf("obisplit -p MICROSAT.PCR_%s_ -t experiment %s_filtered.fastq", 
                            unique(x$library), 
                            unique(x$library))
  system(run.splitlocus)
  # add export/bin to PATH to expose all available tools to the console/shell
  # PATH="$PATH:$HOME/OBITools-1.2.11/export/bin"
  
  # NULL
}, N = N, simplify = FALSE)
