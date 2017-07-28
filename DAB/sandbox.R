library(data.table)

xy.file <- list.files("./DAB/zf_hiseq2/", full.names = TRUE)

xy <- sapply(xy.file, FUN = function(fl) {
  x <- fread(fl)
  samplecols <- names(x)[grepl("sample:", names(x))]
  x[, count := Reduce("+", .SD), .SDcols = samplecols] # make sure correct count
  
  # remove zero alleles
  x <- x[count != 0, ]
  x$count <- NULL
  x$id <- NULL
  x <- melt(x, id.vars = c("seq_length", "sequence"), value.name = "count",
            variable.name = "Sample_Name")
  
  # add information of position, run, and library
  # also, make sample name pretty
  x[, position := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\2", Sample_Name)]
  x[, run := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\3", Sample_Name)]
  x[, Sample_Name := gsub("^sample:(.*)_(\\d+)_PP(\\d+)$", "\\1", Sample_Name)]
  x[, library := gsub("^.*_(DAB\\d+)_.*$", "\\1", basename(fl))]
  x <- setkey(x, run, count)
  x <- x[i = order(-run, count, decreasing = TRUE), j = .(Sample_Name, run, count, seq_length, position, library, sequence), by = .(run)]
  x
}, simplify = FALSE)
xy <- rbindlist(xy)

# find sex based on ZFX/Y
xy[count > 100, j = .N, by = .(Sample_Name, run)]
# N == 2: male
# N == 1: female

i <- xy[, unique(Sample_Name)]
xy[Sample_Name == i[22], .(Sample_Name, run, count, seq_length, sequence)]
xy[Sample_Name == "M05MA", ]
print(xy[Sample_Name == "TM001C", ], 400)
print(xy[sequence == "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccgtcagagctcaaaaagcacatgcg" |
           sequence == "ttgaatcgccaccttttggcggtccacagcaagaactttcctcatatttgtgtggagtgcggtaaaggttttcgtcacccatcagagctcaaaaagcacatgcg", ], 100)
