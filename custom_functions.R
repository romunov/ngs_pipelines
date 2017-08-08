#' For each library will create 8 plates per page. Each plate has 8x12 positions with corresponding sample 
#' name and read count.

showSumsByLibrary <- function(mc) {
  # change function sum to mean, if you're interested in mean
  mcp <- mc[, .(mean.count = sum(Read_Count, na.rm = TRUE)), by = .(Sample_Name, Run_Name, Plate, Position)]
  
  out <- sapply(split(mcp, f = list(mcp$Run_Name)), FUN = function(x) {
    out <- sapply(split(x, f = x$Plate), FUN = function(y) {
      y$xpos <- as.numeric(y$Position)
      xout <- rep(NA, 96)
      lbs <- xout # vector used to fill in labels
      
      xout[y$xpos] <- y[, mean.count]
      lbs[y$xpos] <- sprintf("%s\n%s", y[, Sample_Name], round(y[, mean.count], 0))
      lbs <- matrix(lbs, nrow = 8)
      
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
      
      ggplot(out, aes(x = Var2, y = Var1)) +
        theme_bw() +
        scale_x_discrete(position = "top") +
        theme(axis.ticks = element_blank(), legend.position = "none") +
        xlab("") + ylab("") +
        ggtitle(unique(y$Plate)) +
        geom_tile(aes(fill = value)) +
        geom_text(data = lbs, aes(x = Var2, y = Var1, label = lbl), color = "white")
      
    }, simplify = FALSE)
    
    # plot to file
    pdf(file = sprintf("platecount_%s.pdf", unique(x$Run_Name)), width = 20, height = 20)
    do.call(grid.arrange, c(grobs = out, nrow = 4, as.table = FALSE))
    dev.off()
  }, simplify = FALSE)
  return(invisible(out))
}
