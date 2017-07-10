library(tidyr)
library(uuid)
library(parallel)
library(ggplot2)
library(gridExtra)

minNcount <- 0 # only alleles with  > minNcount will be processed

# This script reads summarized "serie.tab" files and produces genotypes. Certain parts of it are
# run in parallel using the \code{parallel} package which comes shipped with R and works on Unix
# and Windows machines.
# 
# Input files should look like (from microsatTabToseries.py):
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
#
#
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

# On windows system, utilize 4 cores and 46 in Unix. Adapt according to your hardware capacity.
if (Sys.info()["sysname"] == "Windows") {
  ncores <- 4
} else {
  ncores <- 46
}

# Boilerplate stuff to initialize processes.
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(tidyr))

# Find files and exclude possible bad samples (e.g. if you didn't differentiate between POS+1 and POS+2).
xy <- list.files("./data", pattern = "^MICROSAT.*_serie.tab", full.names = TRUE)
xy <- xy[!grepl("_PCR[+-]_", basename(xy))]

message(sprintf("Processing %d files.", length(xy)))

# For each file, read in the data, sort it according to total count of reads, remove unnecessary columns,
# reflow ata into a long format for all repeats, modify or create some new columns like sample, locus, run,
# and finish by reordering the result according to count within each run. This means that first e.g. 6 rows
# will hold six loci for run 1 (P1) ordered by the number of reads.
genotypes <- parSapply(cl = cl, X = xy, FUN = function(x, Nalleles = 6) {
  # genotypes <- sapply(X = xy[grepl("M0P0K.MM", xy)][5], FUN = function(x, Nalleles = 6) {
  # Import data.
  io <- read.table(x, header = TRUE)
  
  # Subset six most common alleles (this can be varied).
  io <- io[order(io$count, decreasing = TRUE), ][1:Nalleles, ]
  
  # Remove unnecessary columns.
  io$id <- NULL # uuid to be added later
  io$count <- NULL # sum of count_run ~ seq_length + series
  
  # Reflow data into a long format.
  io <- gather(io, key = run, value = count_run, 
               c(-sequence, -seq_length, -series))
  
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
  
  # Sort by run and then decreasing count_run.
  io <- io[order(rev(io$run), io$count_run, decreasing = TRUE), ]
  
  io
}, simplify = FALSE)

## Save data into a .RData file in case shit hits the fan.
# save(genotypes, file = "./data/genotypes_gatc_dinalpbear.RData")
# genotypes <- do.call(rbind, genotypes)
# rownames(genotypes) <- NULL
# save(genotypes, file = "./data/genotypes_gatc_dinalpbear_one_df.RData")

load("./data/genotypes_gatc_dinalpbear_one_df.RData")

count.bf <- ggplot(genotypes, aes(x = count_run)) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 500)) +
  ggtitle(sprintf("Read count per sequence", minNcount)) +
  geom_histogram(binwidth = 1)

# Reads of length 0 are nonsense and are removed.
genotypes <- droplevels(genotypes[genotypes$count_run > minNcount, ])

# Save POS data into its own data.frame for tidy purposes.
blkpos <- genotypes[grepl("POS", genotypes$sample), ]
rm(blkpos)
genotypes <- genotypes[!grepl("POS", genotypes$sample), ]
genotypes$sequence <- as.character(genotypes$sequence)
gt <- genotypes # I feel typing "gt" is faster than "genotypes"
rm(genotypes)

# check if each sequence has only one locus
unique(aggregate(locus ~ sequence, data = gt, FUN = function(x) length(unique(x)))$locus)
# ... expecting [1] 1

# Standardize family by locus and length
gt$length <- gsub("^(\\d+)_\\d+", "\\1", gt$allele)

gt <- do.call(rbind, by(data = gt, INDICES = list(gt$locus, gt$length), FUN = function(x) {
  # if (unique(x$locus) == "03" & unique(x$length == "93")) {
  #   browser()
  # }
  # if (unique(x$locus) == "03") browser()
  
  # If allele is detected only once, discard it.
  tmp <- split(x, f = x$sequence)
  num.occur.seq <- sapply(tmp, nrow)
  cand.allel <- names(num.occur.seq[num.occur.seq > 1]) # candidate alleles which occur at least twice
  if (length(cand.allel) == 0) return(NULL)
  x <- x[x$sequence %in% cand.allel, ]
  
  # Sort and re-assign new allele numbers.
  x <- x[order(x$count_run, decreasing = TRUE), ]
  x$new_allele <- as.numeric(factor(x$sequence, levels = unique(x$sequence), labels = 1:length(unique(x$sequence))))
  x$new_allele <- paste(gsub("^(\\d+)_\\d$", "\\1", x$allele), x$new_allele, sep = "_")
  x[order(x$sample, x$run, rev(x$count_run)), ]
})
)
rownames(gt) <- NULL

# order so that sequence is at the end for easier viewing
gt <- gt[, c("sample", "run", "count_run", "locus", "allele", "new_allele", "sequence")]
gt <- gt[order(gt$sample, gt$locus, gt$run, rev(gt$count_run)), ]

# Add run position and tags so that we're able to reconstruct NGS filter.
ngs <- read.table("UA_gatc_8_plate.ngsfilter")
ngs$locus <- gsub("^.*_\\w+_(\\d+|(ZF)+)$", "\\1", ngs[, "V1"])
ngs$position <- gsub("^.*_(\\d+)_P\\d+$", "\\1", ngs$V2)
ngs$sample <- gsub("^(.*)_\\d+_.*$", "\\1", ngs$V2)
ngs$run <- gsub("^.*_\\d+_(P\\d+)$", "\\1", ngs$V2)

ss <- c("sample", "position", "locus", "V3", "run")
ngs <- ngs[!duplicated(ngs[, ss]),][, ss]
names(ngs)[4] <- "tagcombo"

outgt <- merge(x = gt, y = ngs)
outgt$Sample_ID <- UUIDgenerate()

input <- list.files("./rawdata", pattern = "fastq.gz")[1]
outgt$Run_Name <- gsub("_\\d\\.fastq\\.gz", "", input)

outgt$run <- gsub("P", "", outgt$run)

names(outgt) <- c("Sample_Name", "Plate", "Marker", "Read_Count", "old_Allele", "Allele",
                  "Sequence", "Position", "TagCombo", "Sample_ID", "Run_Name")

count.af <- ggplot(outgt, aes(x = Read_Count)) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 500)) +
  ggtitle(sprintf("Read count per sequence after removing \nthose with %s reads or less", minNcount)) +
  geom_histogram(binwidth = 1)

pdf(sprintf("histogram_allele_read_count_%s.pdf", unique(outgt$Run_Name)), width = 4, height = 6)
grid.arrange(count.bf, count.af, ncol = 1)
dev.off()

write.table(outgt, file = "genotypes_dinalpbear_gatc_notrush.txt", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# Find number of different alleles per locus.
apl <- aggregate(Read_Count ~ Marker + Allele, FUN = sum, data = outgt)
apl <- apl[order(apl$Marker, apl$Allele, apl$Read_Count, decreasing = TRUE), ]
write.table(apl, file = "num_alleles_per_locus.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# outgt[gsub("^(\\d+)_(\\d+)$", "\\2", outgt$Allele) == 13, ]
# tmp <- outgt
# tmp$length <- gsub("^(\\d+)_(\\d+)$", "\\1", tmp$Allele)
# tmp$aver <- gsub("^(\\d+)_(\\d+)$", "\\2", tmp$Allele)
# tmp <- tmp[(tmp$length == "98" & tmp$Marker == "14"), ]
# tmp <- split(tmp, f = tmp$Sequence)
# sapply(tmp, write.table, row.names = FALSE, file = "test.txt", append = TRUE)
# 
# unique(tmp[tmp$length == "98"& tmp$Marker == "68", "Sequence"])
# tmp[tmp$length == "98"& tmp$Marker == "68", "Allele"]

