# This script will output raw genotypes for ZF.
# Input is the folder where _ZF.uniq.tab files (as produced by the last command in the pipeline - obitab).

library(data.table)

#### INPUT ####
# point this line to the folder where .uniqe.tab files of ZF are located
xy.file <- list.files("./DAB/zf_hiseq2/", full.names = TRUE)
# name of the file into which the results are to be written
file.out <- "./DAB/data/genetic_sex_hiseq2.txt"
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
sexy[, sex := ""]
sexy[sequence %in% seqUAX, sex := "X"]
sexy[sequence %in% seqUAY, sex := "Y"]

# write intermediate (raw) result
fwrite(sexy, file = file.out)

x <- sexy[Sample_Name %in% "M00CL", -"sequence"]
x <- x[-2]
x <- x[order(Sample_Name, run, sex)]
x[, xrle := seq(1:.N), by = .(run)]
x <- dcast(x, run + Sample_Name + seq_length + position + library ~ xrle, value.var = c("count", "sex"))
x[, ratio := count_1/count_2]
x
