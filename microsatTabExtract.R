microsatTabExtract <- function(filename, sampleName) {
  #### extract data according to sample name and print these data ordered by counts
  ## the script keeps the two first columns : id and count (updated by the script)
  ## and the last two columns : seq_length and sequence
  ## all other columns that do not match with the sample name prefix provided are removed in the output
  
  #### collect and check arguments provided
  # args <- commandArgs(TRUE)
  # 
  # if (length(args) < 2 | length(args) > 2)
  # {
  #   print ("##### Error Message #####")
  #   print ("you need to provide two arguments in this order :")
  #   print ("1- the name of the input tab file")
  #   print ("2- the abbreviated name used in the file for the samples (ex : CLF or UAS")
  #   stop()
  # }
  
  # filename <- args[1]
  # sampleName <- args[2]
  
  #### read file and get data
  data <- read.table(filename, header = TRUE)
  nom <- colnames(data)
  ext <- sapply(nom, function(x) strsplit(x,'_')[[1]][1])
  
  extSamplesNames <- sapply(ext, function(x) substr(x, 8, nchar(x)))
  sampleColNumbers <- which(extSamplesNames == sampleName)
  id <- data$id
  sequence <- data$sequence
  seq_length <- data$seq_length
  dfExt <- data[, sampleColNumbers]
  dfExtColNames <- colnames(dfExt)
  
  ext2 <- as.numeric(sapply(dfExtColNames, function(x) strsplit(x, '_')[[1]][2]))
  sampleColNumbers2 <- which(ext2 > 0)
  dfExt2 <- dfExt[, sampleColNumbers2]
  count <- apply(dfExt2, MARGIN = 1, FUN = sum)
  
  #### create dataframe and print to output file
  # df <- as.data.frame(cbind(id, count, dfExt2, seq_length, sequence))
  df <- data.frame(id, count, dfExt2, seq_length, sequence)
  df <- df[order(df$count, decreasing = TRUE), ]
  fileoutname <- paste("MICROSAT.PCR", "_", sampleName, "_", strsplit(filename, "_")[[1]][4], sep = "")
  
  write.table(df, fileoutname, quote = FALSE, row.names = FALSE, sep = "\t")
}
