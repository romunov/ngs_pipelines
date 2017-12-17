# This script will extract alleles X and Y and call genetic sex. It produces 2-3 files, depending on your wishes.
# Go through comments to see the logic of the way genetic sex is called.
# Input is the folder where _ZF.uniq.tab files (as produced by the last command in the pipeline - obitab).

library(data.table)
source("custom_functions.R")

#### INPUT ####
# Point this line to the folder where .uniqe.tab files of ZF are located.
# xy.file <- list.files("./DAB/zf_hiseq1/", full.names = TRUE)
# ngs.files <- list.files("./DAB/1_ngsfilters_hiseq1/", pattern = ".ngsfilter", full.names = TRUE)
xy.file <- list.files("./DAB_GATC2/zf", full.names = TRUE)
ngs.files <- list.files("./DAB_GATC2/1_ngsfilters/", pattern = ".ngsfilter", full.names = TRUE)

# Name of the file(s) into which the results are to be written. There is one more 
# possible intermediate file towards the end.
# Switch between HISEQ1 and HiSEQ2 runs.
# file.called <- "./DAB/data/dab_called_genetic_sex_hiseq1.txt"
# file.freqs <- "./DAB/data/frequency_of_sequences_by_marker_hiseq1_sex.txt"

file.called <- "./DAB_GATC2/data/dab_called_genetic_sex_gatc2.txt"
file.freqs <- "./DAB_GATC2/data/frequency_of_sequences_by_marker_gatc2_sex.txt"
countTS <- 20 # reads that do not match exact sequence should have this number of reads,
              # otherwise they are discarded
lowCount <- 100 # if unrecognized as sex sequence and below this threshold, flag as "L"
junkTS <- 0.15 # when scanning for junk, if a sequence is less than this of max sequence, discard
dbTS <- 0.5 # if sequence is below dbTS relative to max sequence, flag as disbalanced
#### INPUT ####

# Exact sequences from X and Y chromosomes
seqUAX <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccgtcagagctcaaaaagcacatgcg"
seqUAY <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccatcagagctcaaaaagcacatgcg"

# For each file/sample, extract alleles along with Position, run, sample name, Run_Name and order properly.
system.time(xy <- sapply(xy.file, FUN = function(fl) {
  out.cols <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", 
                "length", "Position", "called", "flag", "stutter", "Sequence")
  # read in the data
  x <- fread(fl)
  
  # find sample names
  samplecols <- names(x)[grepl("sample:", names(x))]
  # make sure total count is correct sum of all samples
  x[, count := Reduce("+", .SD), .SDcols = samplecols] 
  
  # remove zero alleles and unneeded columns
  x <- x[count != 0, ]
  x[, c("count", "id") := NULL]
  names(x)[grepl("sequence", names(x))] <- "Sequence"
  names(x)[grepl("seq_length", names(x))] <- "length"
  
  # reshape the data from wide to long format
  x <- melt(x, id.vars = c("length", "Sequence"), value.name = "Read_Count",
            variable.name = "Sample_Name")
  
  # add information all extra information
  x[, Position := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\2", Sample_Name)]
  x[, Plate := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\3", Sample_Name)]
  x[, Sample_Name := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\1", Sample_Name)]
  # x[, Run_Name := gsub("^.*_(DAB\\d+)_.*$", "\\1", basename(fl))]
  x[, Run_Name := sprintf("DAB%s", gsub("^MICROSAT\\.PCR-(\\d+)_.*$", "\\1", basename(fl)))]
  x[, Marker := "ZF"]
  x[, called := TRUE]
  x[, flag := ""]
  x[, stutter := FALSE]
  x[, TagCombo := ""]

  x <- x[, out.cols, with = FALSE]
  
  # put everything in proper order
  x <- x[i = order(-Plate, Read_Count, decreasing = TRUE), 
         j = out.cols, with = FALSE, 
         by = .(Plate)]
  x
}, simplify = FALSE))

# Combine the result into one data.table and exclude zero counts.
xy <- rbindlist(xy)
xy <- xy[Read_Count > 0, ]
xy <- xy[length >= 100, ] # retain only sequences longer or equal to 100

# Sequences of X and Y, sent in by Stéphane Lobreaux (Université Grenoble Alpes). Use them to replace with
# human readable format.
uaxy <- c(seqUAX, seqUAY)

# Extract only sequences with sufficient count.
sexy <- xy[Read_Count > countTS, ]

#### add TagCombo ####
tc <- sapply(ngs.files, FUN = read.table, simplify = FALSE)
tc <- rbindlist(tc)
rownames(tc) <- NULL
tc <- as.data.table(tc[, c(1, 2, 3)])
names(tc) <- c("V1", "V2", "TagCombo")

# create columns by which to merge
sexy[, fn := sprintf("%s_%s_PP%s", Sample_Name, Position, Plate)]
sexy[, fl := sprintf(sprintf("UA_MxRout1_%s", Marker))]

sexy <- merge(sexy, tc, by.x = c("fn", "fl"), by.y = c("V2", "V1"))
sexy[, fn := NULL]
sexy[, fl := NULL]
#####
# keep only sequences with +-5 of length of X
sexy <- sexy[length < nchar(seqUAX) + 5 & length > nchar(seqUAX) - 5, ]

# if sequence is not known to be from UAX or UAY, mark as low if below
# a given threshold
sexy[!(Sequence %in% uaxy) & Read_Count <  lowCount, flag := paste(flag, "L", sep = "")]

# clean based on function cleanZF
sps <- split(sexy, f = list(sexy$Sample_Name, sexy$Plate))
system.time(spsclean <- sapply(sps, FUN = cleanZF, ts = junkTS, db = dbTS, simplify = FALSE))
spsclean <- rbindlist(spsclean)
spsclean[Sample_Name %in% sample(spsclean$Sample_Name, 1), -c("Sequence", "TagCombo")]

fwrite(spsclean[order(Sample_Name, Plate, Read_Count)], file = file.called, sep = "\t")

freq.of.seq.by.marker <- spsclean[, .N, by = .(Marker, Sequence)][order(-Marker, N, decreasing = TRUE), .(Marker, N, Sequence)]
fwrite(freq.of.seq.by.marker, file = file.freqs, sep = "\t")

# number of flagged samples:
spsclean[, .(numsamples = length(unique(Sample_Name))), by = .(flag)]
