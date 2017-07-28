#' Extract data for a given sample from "locus reservoir" file.
#' 
#' @param filename Character. Path to the filename which holds data on all samples for individual locus.
#' @param samplename Character. Sample name.
#' @param outdir Path where the result is to be stored. No trailing slash. Defaults to "." (here).
#' 
microsatTabExtract <- function(filename, samplename, outdir = ".") {
  xy <- fread(filename, stringsAsFactors = FALSE)
  
  # prepare library designation
  lib <- gsub("^MICROSAT.*_JFV-(\\d+)_.*$", "\\1", basename(filename))
  lib <- sprintf("DAB%02d", as.numeric(lib))
  
  fileoutname <- sprintf("%s/MICROSAT.PCR_%s_%s_%s.uniq.tab",
                         outdir,
                         lib,
                         samplename,
                         # extract locus name
                         gsub("^MICROSAT\\.PCR_.*_([[:alnum:]]+)\\.uniq\\.tab$", "\\1", basename(filename)))

  # From a file which holds data for all samples and all replications, extract only the desired file.
  samcols <- colnames(xy)[grepl(samplename, colnames(xy))]
  
  if (length(samcols) == 0) {
    message(sprintf("Found no alleles for %s", fileoutname))
    return(NULL)
  }
  subcols <- c("id", "count", samcols, "seq_length", "sequence")
  # Next to sample data, extract also id, count of reads, length of the sequence and the sequence itself.
  out <- xy[, ..subcols]
  
  # Sort by count in decreasing order.
  out$count <- rowSums(out[, samcols, with = FALSE]) # recalculate count for this sample
  out <- out[order(out$count, decreasing = TRUE), ]
  
  write.table(out, file = fileoutname, quote = FALSE, row.names = FALSE, sep = "\t")
}
