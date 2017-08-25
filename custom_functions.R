#' For each library will create 8 plates per page. Each plate has 8x12 positions with corresponding sample 
#' name and read count of all sequences.
#' Use parameter `loc` to output result to a specific folder. Do not forget a trailing slash. Defaults to `getwd()`.
#' 
#' Parameter pos will direct how plates are position in the output. `by_plate` it will preserve the order where
#' plates are 1-4 in the first column and 5-8 in the second column. `by_sample` will organize plates to have
#' plates with same replication together. `pos` can also be a numeric integer vector of the same length as there
#' are number of plates. This will designated in which order the plates are to be plotted. Keep in mind that 
#' the plotting order is column-wise.
#' 
showSumsByLibrary <- function(mc, loc = "./", pos = c("by_plate", "by_sample")) {
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  # change function sum to mean, if you're interested in mean
  mcp <- mc[, .(mean.count = sum(Read_Count, na.rm = TRUE)), by = .(Sample_Name, Run_Name, Plate, Position)]
  
  out <- sapply(split(mcp, f = list(mcp$Run_Name)), FUN = function(x) {
    out <- sapply(split(x, f = x$Plate), FUN = function(y) {
      # create a vector and add data for correct positions
      # this is necessary to preserve missing positions
      y$xpos <- as.numeric(y$Position)
      xout <- rep(NA, 96)
      lbs <- xout # vector used to fill in labels
      nms <- xout # vector used to fill in names
      
      xout[y$xpos] <- y[, mean.count]
      lbs[y$xpos] <- sprintf("%s\n%s", y[, Sample_Name], round(y[, mean.count], 0))
      lbs <- matrix(lbs, nrow = 8)
      
      nms[y$xpos] <- y[, Sample_Name]
      nms <- matrix(nms, nrow = 8)
      colnames(nms) <- 1:12
      rownames(nms) <- LETTERS[1:8]
      nms <- reshape2::melt(nms)
      nms$isnegcont <- NA
      sn <- y[grepl("(^AC\\..*$|^AK\\..*$|NeKo.*$)", Sample_Name), Sample_Name]
      nms[nms$value %in% sn, "isnegcont"] <- TRUE
      
      out <- matrix(xout, nrow = 8)
      colnames(out) <- 1:12
      rownames(out) <- LETTERS[1:8]
      out <- reshape2::melt(out)
      out$Var2 <- as.factor(out$Var2)
      out$Var1 <- factor(out$Var1, levels = rev(levels(out$Var1)), ordered = TRUE)
      
      rownames(lbs) <- LETTERS[1:8]
      colnames(lbs) <- 1:12
      lbs <- melt(lbs, value.name = "lbl")
      lbs$Var2 <- as.factor(lbs$Var2)
      lbs$Var1 <- factor(lbs$Var1, levels = rev(levels(lbs$Var1)), ordered = TRUE)
      
      out <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), list(out, nms, lbs))
      
      nk <- sapply(na.omit(out[out$isnegcont == TRUE, c("Var1", "Var2")]), as.numeric, simplify = FALSE)
      nk <- as.data.frame(nk)
      
      ggplot(out, aes(x = Var2, y = Var1)) +
        theme_bw() +
        scale_x_discrete(position = "top") +
        theme(axis.ticks = element_blank(), legend.position = "none") +
        xlab("") + ylab("") +
        ggtitle(unique(y$Plate)) +
        geom_tile(aes(fill = value.x)) +
        geom_rect(data = nk, aes(xmin = Var2 - 0.5, xmax = Var2 + 0.5, ymin = Var1 - 0.5, ymax = Var1 + 0.5), 
                  fill = NA, color = "red", size = 1) +
        geom_text(aes(label = lbl), size = 3.5, color = "white")
      
    }, simplify = FALSE)
    
    # decide if plates will be ordered to check plate performance or sample performance
    if (pos == "by_sample") out <- out[c(1, 7, 3, 5, 2, 8, 4, 6)]
    if (is.numeric(pos)) {
      # stop if pos has insufficient number of positions specified
      stopifnot(length(pos) == length(out))
      out <- out[pos]
    }
    
    # plot to file, loc should have a trailing slash included
    if (!grepl(".*/$", loc)) {
      stop("Trailing slash in `pos` not found. Please check your path.")
    }
    
    pdf(file = sprintf("%splatecount_%s.pdf", loc, unique(x$Run_Name)), width = 20, height = 20)
    do.call(grid.arrange, c(grobs = out, nrow = 4, as.table = FALSE))
    dev.off()
  }, simplify = FALSE)
  return(invisible(out))
}

#' Clean ZF alleles.
#' @param fb Data as organized from the ngs pipeline.
#' @param ts Threshold value under which sequences are discarded relative to the sequence
#' with the highest number of reads in one sample * plate combination.

cleanZF <- function(fb, ts, db) {
  # browser()
  if (nrow(fb) == 0) return(NULL)
  # omit all sequences that do not reach ts% of the highest read.  
  lvls <- fb$Read_Count/max(fb$Read_Count)
  fb <- fb[lvls > ts, ]
 
  # if an allele is significantly below a db threshold, flag it as disbalanced
  lvls <- fb$Read_Count/max(fb$Read_Count)
  fb[lvls < db, ]$flag <- paste(fb[lvls < db, ]$flag, "D", sep = "")
  # fb[lvls < db, flag := paste(flag, "D", sep = "")]
  fb
}
