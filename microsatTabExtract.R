#' Extract data for a given sample from "locus reservoir" file.
#' 
#' @param filename Character. Path to the filename which holds data on all samples for individual locus.
#' @param samplename 
microsatTabExtract <- function(filename, samplename) {
  xy <- read.table(filename, header = TRUE)
  
  # From a file which holds data for all samples and all replications, extract only the desired file.
  samcols <- colnames(xy)[grepl(samplename, names(xy))]
  # Next to sample data, extract also id, count of reads, length of the sequence and the sequence itself.
  out <- xy[, c("id", "count", samcols, "seq_length", "sequence")]
  fileoutname <- sprintf("MICROSAT.PCR_%s_%s.uniq.tab",
                         samplename,
                         # extract locus name
                         gsub("^MICROSAT\\.PCR_.*_([[:alnum:]]+)\\.uniq\\.tab$", "\\1", basename(filename)))
  
  # Sort by count in decreasing order.
  out <- out[order(out$count, decreasing = TRUE), ]
  out$count <- rowSums(out[, samcols]) # recalculate count for this sample
  
  write.table(out, file = fileoutname, quote = FALSE, row.names = FALSE, sep = "\t")
}
