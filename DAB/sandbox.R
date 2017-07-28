library(data.table)
library(tidyr)

xy.file <- "./DAB/zf_hiseq2/MICROSAT.PCR_DAB13_M00CL_ZF.uniq.tab"
xy <- fread(xy.file,
                 header = TRUE)

samplecols <- names(xy)[grepl("sample:", names(xy))]

# make sure the sum is correct
xy[, count := Reduce("+", .SD), .SDcols = samplecols]

# remove zero alleles
xy <- xy[count > 0, ]
xy$count <- NULL
xy$id <- NULL

xy <- melt(xy, id.vars = c("seq_length", "sequence"), value.name = "count",
          variable.name = "Sample_Name")

# add information of position, run, and library
# also, make sample name pretty
xy[, position := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\2", Sample_Name)]
xy[, run := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\3", Sample_Name)]
xy[, Sample_Name := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\1", Sample_Name)]
xy[, library := gsub("^.*_(DAB\\d+)_.*$", "\\1", basename(xy.file))]


     
