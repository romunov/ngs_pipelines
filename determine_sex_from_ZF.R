# This script will extract alleles X and Y and call genetic sex. It produces 2-3 files, depending on your wishes.
# Go through comments to see the logic of the way genetic sex is called.
# Input is the folder where _ZF.uniq.tab files (as produced by the last command in the pipeline - obitab).

library(data.table)

#### INPUT ####
# Point this line to the folder where .uniqe.tab files of ZF are located.
xy.file <- list.files("./DAB/zf_hiseq2/", full.names = TRUE)
# xy.file <- list.files("./DAB/zf_hiseq1/", full.names = TRUE)
# Name of the file(s) into which the results are to be written. There is one more 
# possible intermediate file towards the end.
# Switch between HISEQ1 and HiSEQ2 runs.
file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq1.txt"
file.called <- "./DAB/data/dab_called_genetic_sex_hiseq1.txt"
# file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq2.txt"
# file.called <- "./DAB/data/dab_called_genetic_sex_hiseq2.txt"
#### INPUT ####

# For each file/sample, extract alleles along with position, run, sample name, library and order properly.
xy <- sapply(xy.file, FUN = function(fl) {
  # read in the data
  x <- fread(fl)
  
  # find sample names
  samplecols <- names(x)[grepl("sample:", names(x))]
  # make sure total count is correct sum of all samples
  x[, count := Reduce("+", .SD), .SDcols = samplecols] 
  
  # remove zero alleles and unneeded columns
  x <- x[count != 0, ]
  x[, c("count", "id") := NULL]
  
  # reshape the data from wide to long format
  x <- melt(x, id.vars = c("seq_length", "sequence"), value.name = "count",
            variable.name = "Sample_Name")
  
  # add information of position, run, and library
  # also, make sample name pretty
  x[, position := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\2", Sample_Name)]
  x[, run := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\3", Sample_Name)]
  x[, Sample_Name := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\1", Sample_Name)]
  x[, library := gsub("^.*_(DAB\\d+)_.*$", "\\1", basename(fl))]

  # put everything in proper order
  x <- x[i = order(-run, count, decreasing = TRUE), 
         j = .(Sample_Name, count, seq_length, position, library, sequence), 
         by = .(run)]
  x
}, simplify = FALSE)
# Combine the result into one data.table and exclude zero counts.
xy <- rbindlist(xy)
xy <- xy[count > 0, ]

# Sequences of X and Y, sent in by Stéphane Lobreaux (Université Grenoble Alpes). Use them to replace with
# human readable format.
seqUAX <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccgtcagagctcaaaaagcacatgcg"
seqUAY <- "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccatcagagctcaaaaagcacatgcg"
uaxy <- c(seqUAX, seqUAY)

# Extract only 100% sequences for X and Y.
sexy <- xy[sequence %in% uaxy, ]

# Add human readable designations for X and Y.
sexy[, seq := ""]
sexy[sequence %in% seqUAX, seq := "X"]
sexy[sequence %in% seqUAY, seq := "Y"]

# write intermediate (raw) result
fwrite(sexy, file = file.out.seq)

# Raw sequence is of no importance (from here on), so we remove it. It can always be
# reconstructed from X and Y
sexy[, sequence := NULL]

# If you are skipping the above steps, you can always read in file.out.seq and go from here.
# sexy <- fread(file.out.seq)

# Find samples with more than two alleles which may be causing problems. Uncommend fwrite line to
# output this to a readable file. Some samples may be salveageable, but only manually.
tm <- sexy[Sample_Name %in% sexy[, .N, by = .(Sample_Name, run, library)][N > 2, Sample_Name], ]
tm <- tm[order(Sample_Name, run, library), ]
# fwrite(tm, file = "./DAB/data/samples_with_more_than_two_alleles.txt", sep = "\t")

# Exclude samples with potential flaws - handle them manually if critical.
sexy <- sexy[!(Sample_Name %in% unique(tm$Sample_Name)), ]

# Recast the data so that X, Y and count go from long to wide format.
# before:        after:
# seq count         
#   X   513       X  Y
#   Y     1     513  1
# It is fairly easy to read genotypes in this format.
sexy <- dcast(sexy, run + Sample_Name + seq_length + position + library ~ seq, 
      value.var = c("count"))

# Since some may not have Y (or X), we turn NAs to 0. This and ordering may somewhat
# improve readability.
sexy <- sexy[order(library, Sample_Name, run)]
sexy[, X := ifelse(is.na(X), 0, X)]
sexy[, Y := ifelse(is.na(Y), 0, Y)]

# It would appear males have either balanced X and Y or skewed number of reads for
# Y. Females have very skewed number of reads toward X. We use this piece of
# information to call genetic sex. Using coefficient of variation we judge whether
# X or Y is significantly higher than the other (dispersed far from the mean). If 
# not, it indicates balanced X and Y and not otherwise.

# Calculate coefficient of variation - this will tell us if the difference
# between X and Y is "large" or not. If you plot cv, you will notice that it 
# clusters around zero and above 1. This can be used to define a cut-off.
sexy[, mean := apply(sexy[, .(X, Y)], 1, mean)]
sexy[, cv := apply(sexy[, .(X, Y)], 1, sd) / mean]

# Small values indicate small differences between X and Y --> calling it male.
cutoff.cv <- 0.75
# Define a threshold below which absolute number of reads is considered as 
# untrustworthy and thus forloring genetic sex calling.
cutoff.abs <- 5

# Deem male:
# If coefficient of variation is very small, this indicates X and Y are balanced.
sexy[cv < cutoff.cv, sex := "M"] # this also takes care of cases where X == Y
# If Y is amplified more than X, deem this male.
sexy[Y > mean & X < mean & cv > cutoff.cv, sex := "M"]
# Deem female:
# If X is significantly higher than Y, deem this female.
sexy[X > mean & Y < mean & cv > cutoff.cv, sex := "F"]
# ID of samples with really weak, reads are deemed unreliable.
sexy[X < cutoff.abs & Y < cutoff.abs, sex := NA]

# This just rounds CV to improve readability.
sexy[, cv := round(cv, 2)]
# Order before finally writing the result to disk.
sexy <- sexy[, .(Sample_Name, run, position, library, mean, cv, X, Y, sex)]
fwrite(sexy, file = file.called, sep = "\t")
