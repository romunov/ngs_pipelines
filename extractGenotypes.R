#' This function extracts genotypes into a uniform form.
extractGenotypes <- function(x) {
  rxstr <- "(^ctable2n_MICROSAT\\.PCR)_(.*)_(\\d+)_(serie\\.tab$)"
  samplename <- gsub(rxstr, "\\2", x)
  locusname <- gsub(rxstr, "\\3", x)
  
  y <- read.table(x, header = TRUE)
  find.nas <- apply(y, MARGIN = 1, function(x) all(is.na(x) | x == 0))
  rn <- rownames(y)
  
  out <- y[!find.nas,]
  colorder <- c("samplename", "locusname", "run", paste("allele", 1:6, sep = ""),
                paste("count_allele", 1:6, sep = ""))
  
  if (nrow(out) == 0) {
    ot <- data.frame(matrix(rep(NA, ncol(out)), nrow = 1))
    colnames(ot) <- colnames(out)
    out <- rbind(out, ot)
    out$samplename <- samplename
    out$locusname <- locusname
    rownames(out) <- out$samplename
    out$run <- NA
    return(out[, colorder])
  }
  
  rn <- rownames(out)
  rnstr <- "(^sample)\\.(.*)_(\\d+)_(.*)$"
  out$run <- gsub(rnstr, "\\4", rn)
  out$samplename <- samplename
  out$locusname <- locusname
  
  # remove consensus genotypes from the output
  out <- out[!(out$run %in% c("CNS")), ]
  rownames(out) <- NULL
  out[, colorder]
}
