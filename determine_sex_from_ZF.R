# This script will output raw genotypes for ZF.
# Input is the folder where _ZF.uniq.tab files (as produced by the last command in the pipeline - obitab).

library(data.table)
library(ggplot2)
library(plotly)

#### INPUT ####
# point this line to the folder where .uniqe.tab files of ZF are located
xy.file <- list.files("./DAB/zf_hiseq2/", full.names = TRUE)
# xy.file <- list.files("./DAB/zf_hiseq1/", full.names = TRUE)
# name of the file into which the results are to be written
# file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq1.txt"
# file.called <- "./DAB/data/dab_called_genetic_sex_hiseq1.txt"
file.out.seq <- "./DAB/data/dab_genetic_sex_with_sequence_hiseq2.txt"
file.called <- "./DAB/data/dab_called_genetic_sex_hiseq2.txt"
#### INPUT ####

# for each file, extract alleles, along with position, run, sample name, library and order properly
xy <- sapply(xy.file, FUN = function(fl) {
  x <- fread(fl)
  samplecols <- names(x)[grepl("sample:", names(x))]
  x[, count := Reduce("+", .SD), .SDcols = samplecols] # make sure count is correct
  
  # remove zero alleles and unneeded columns
  x <- x[count != 0, ]
  x[, c("count", "id") := NULL]
  
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
# combine the result into one data.table
xy <- rbindlist(xy)
# exclude zero counts
xy <- xy[count > 0, ]

# sequences of X and Y, sent in by Stephane L.
seqUAX="ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccgtcagagctcaaaaagcacatgcg"
seqUAY="ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccatcagagctcaaaaagcacatgcg"
uaxy <- c(seqUAX, seqUAY)

# extract only correct sequences for X and Y
sexy <- xy[sequence %in% uaxy, ]

# add human readable designations for X and Y
sexy[, seq := ""]
sexy[sequence %in% seqUAX, seq := "X"]
sexy[sequence %in% seqUAY, seq := "Y"]

# write intermediate (raw) result
fwrite(sexy, file = file.out.seq)

sexy[, sequence := NULL]

# sexy <- fread(file.out.seq)

# samples with more than two alleles
tm <- sexy[Sample_Name %in% sexy[, .N, by = .(Sample_Name, run, library)][N > 2, Sample_Name], ]
tm <- tm[order(Sample_Name, run, library), ]
# fwrite(tm, file = "./DAB/data/samples_with_more_than_two_alleles.txt", sep = "\t")

# exclude samples with potential flaws
sexy <- sexy[!(Sample_Name %in% unique(tm$Sample_Name)), ]

# find which samples have more than two alleles per run
sexy <- dcast(sexy, run + Sample_Name + seq_length + position + library ~ seq, 
      value.var = c("count"))

sexy <- sexy[order(library, Sample_Name, run)]
sexy[, X := ifelse(is.na(X), 0, X)]
sexy[, Y := ifelse(is.na(Y), 0, Y)]

# calculate coefficient of variation - this will tell us if the difference
# between X and Y is large or not
sexy[, mean := apply(sexy[, .(X, Y)], 1, mean)]
sexy[, cv := apply(sexy[, .(X, Y)], 1, sd) / mean]

cutoff.cv <- 0.75 # small means difference between X and Y are small and vice versa
cutoff.abs <- 5 # how many reads should there be before we reliably determine sex

# if Y is amplified more than X, deem this male
sexy[cv < cutoff.cv, sex := "M"] # this also takes care of cases where X == Y
sexy[Y > mean & X < mean & cv > cutoff.cv, sex := "M"]
# if X is significantly higher than Y, deem this female
sexy[X > mean & Y < mean & cv > cutoff.cv, sex := "F"]
# ID of samples with really weak reads are deemed unreliable
sexy[X < cutoff.abs & Y < cutoff.abs, sex := NA]

sexy[, cv := round(cv, 2)]
sexy <- sexy[, .(Sample_Name, run, position, library, mean, cv, X, Y, sex)]
fwrite(sexy, file = file.called, sep = "\t")
