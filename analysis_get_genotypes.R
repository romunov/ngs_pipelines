library(tidyr)
library(uuid)
library(parallel)

# This script reads summarized "serie.tab" files and produces genotypes. The script can be
# run in parallel.
# 
# Input should look like (from microsatTabToseries.py):
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
#  And genotype output is:
#       sample run count_run locus allele
# 3  M0P0K.MM  P1         9   064   97_1
# 9  M0P0K.MM  P2        88   064   97_1
# 15 M0P0K.MM  P3         5   064   97_1
# 21 M0P0K.MM  P4        51   064   97_1
# 27 M0P0K.MM  P5        17   064   97_1
# 33 M0P0K.MM  P6        92   064   97_1
# 39 M0P0K.MM  P7        16   064   97_1
# 45 M0P0K.MM  P8        66   064   97_1
# sequence
# 3  gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 9  gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 15 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 21 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 27 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 33 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 39 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
# 45 gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
#'  

if (Sys.info()["sysname"] == "Windows") {
  ncores <- 4
} else {
  ncores <- 46
}

cl <- makeCluster(ncores)
clusterEvalQ(cl, library(tidyr))

xy <- list.files("./data", pattern = "^MICROSAT.*_serie.tab", full.names = TRUE)
xy <- xy[!grepl("_PCR[+-]_", basename(xy))]

message(sprintf("Processing %d files.", length(xy)))

genotypes <- parSapply(cl = cl, X = xy, FUN = function(x, Nalleles = 6) {
# genotypes <- sapply(X = xy, FUN = function(x, Nalleles = 6) {
  # Import data.
  io <- read.table(x, header = TRUE)
  # Subset six most common alleles (this can be varied).
  io <- io[order(io$count, decreasing = TRUE), ][1:Nalleles, ]
  # Reflow data into a long format.
  io <- gather(io, key = run, value = count_run, 
               c(-id, -count, -sequence, -seq_length, -series))
  
  # Remove unnecessary columns.
  io$id <- NULL # uuid to be added later
  io$count <- NULL # sum of count_run ~ seq_length + series
  
  # Add sample names, locus, allele and run number.
  # Explanation of regex:
  # string should begin with "sample.", then find a *group*, then find _,
  # then find digits, followed by _, then by P, followed by more digits
  # and finally the string should end on that. Extract *group*.
  # The next two regex statements are similar, except they vary in number
  # of position of groups. \\d{1} means "find 1 digit".
  io$sample <- gsub("^.*_(.*)_(\\d+)_serie.tab$", "\\1", basename(x))
  io$locus <- gsub("^.*_(.*)_(\\d+)_serie.tab$", "\\2", basename(x))
  io$run <- gsub("(^sample\\..*)_(P\\d{1})", "\\2", io$run)
  io$allele <- paste(io$seq_length, io$series, sep = "_")
  
  # Notice that `serie` and `seq_length` are no longer included as they are
  # encapsulated in `allele`.
  io <- io[, c("sample", "run", "count_run", "locus", "allele", "sequence")]
  
  # # TODO: This sanity check should be omitted.
  # # Sanity check that each allele contains only one kind of a sequence.
  # sanch <- split(io, f = io$allele)
  # sanch <- sapply(sanch, FUN = function(sc) {
  #   # If per allele there is only one sequence, return TRUE, otherwise FALSE.
  #   if (length(unique(sc$sequence)) == 1) return(TRUE) else return(FALSE)
  # })
  # if (all(sanch)) {
  #   message(sprintf("Alleles %s locus %s OK.", unique(io$sample), unique(io$locus)))
  # } else {
  #   message("Shit @ %s locus %s.", unique(io$sample), unique(io$locus))
  # } 
  io
}, simplify = FALSE)

# Genotype processing. For instance:
# * removing series which have less than a threshold of reads
