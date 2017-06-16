### script to analyse microsat NGS tab file per sample through sequence family
microsatAlleles <- function(nrep, motifLength, numsamples, numseqfamily, input) {
  nrep <- as.numeric(nrep)
  motifLength <- as.numeric(motifLength)
  NbEch <- as.numeric(numsamples)
  minNbOfSeqPerFamily <- as.numeric(numseqfamily)
  filename <- input
  
  ### read csv data file and remove def column
  data <- read.table(filename, head = TRUE)

  ######### filtering according to the minimal number of sequences per family provided
  familySeqCount <- tapply(data$count, data$series, sum)
  series <- unique(data$series)
  familydf <- as.data.frame(cbind(familySeqCount, series))
  seriesok <- familydf$series[familydf$familySeqCount > minNbOfSeqPerFamily]
  data <- data[is.element(data$series, seriesok), ]
  
  ### def first and last sample positions (sample...R.) (voir regexp)
  # def repet tags to use according to repet number
  Rrep <- paste("P", 1:nrep, sep = "")
  
  # sort data according to seq length values
  data <- data[order(data$seq_length), ]
  
  ### get sample numbers and experiment name
  df <- data[, colnames(data)[grepl("_P\\d+$", colnames(data))]]
  
  nom <- colnames(df)
  dataColNb <- dim(data)[2]
  ext <- as.numeric(sapply(nom, function(x) strsplit(x, '_')[[1]][2]))
  sampleNb <- unique(ext)
  nExt <- strsplit(nom[1], '_')[[1]][1]
  nameExp <- substr(nExt, 8, nchar(nExt))
  
  ###### prepare table2 to transfer genotype data id + 6 geno + 6 counts + CNS
  ### create and initialize table 2
  table2 <- as.data.frame(matrix(data = NA, nrow = NbEch * (nrep + 1), ncol = 12))

  namesTable2 = c()
  for (i in 1:(NbEch)) {
    samp <- i
    add <- ""
    if (samp < 10) {
      add = "00"
    }
    if (samp > 9 & samp < 100) {
      add <- "0"
    }
    
    for (j in 1:(nrep)) {
      currentName <- paste("sample.", nameExp, '_', add, strtoi(samp), '_R', strtoi(j), sep = "")
      namesTable2 <- c(namesTable2, currentName)
    }
    cnsName <- paste("sample.", nameExp, '_', add, strtoi(samp), '_CNS', sep = "")
    namesTable2 <- c(namesTable2, cnsName)
  }
  
  rownames(table2) <- namesTable2
  colNamesTable2 <- c("allele1","allele2","allele3","allele4","allele5","allele6","count_allele1","count_allele2","count_allele3","count_allele4","count_allele5","count_allele6")
  colnames(table2) <- colNamesTable2
  
  ### loop over sample columns
  dfColNames <- colnames(df)
  indice <- 1
  begin <- 1
  end <- 1
  # init sampNumber
  currentSampNumber <- as.numeric(strsplit(dfColNames[1], '_')[[1]][2])
  
  while (indice <= (length(dfColNames) + 1)) {
    if (indice <= (length(dfColNames))) {
      sampNumber <- as.numeric(strsplit(dfColNames[indice], '_')[[1]][2])
    }
    else {
      sampNumber <- "" ###  to analyze the last sample
    }
    if (currentSampNumber != sampNumber & indice != 1) { ### ajouter ici une condition pour traiter le dernier echantillon ou==length()
      ### extract corresponding data
      sampleDf <- df[, begin:end]
      if (begin == end) {
        sampRepet <- 1
        sum <- sampleDf
        serie <- data$serie
        seq_length <- data$seq_length
        sample <- sampleDf
        sampleDf <- as.data.frame(cbind(sample, serie, seq_length))
        namesCol <- c(dfColNames[indice - 1], "serie", "seq_length")
        colnames(sampleDf) <- namesCol
      } else {
        sampRepet <- dim(sampleDf)[2]
        sum <- apply(sampleDf, MARGIN = 1, FUN = sum)
        serie <- data$serie
        seq_length <- data$seq_length
        sampleDf <- as.data.frame(cbind(sampleDf, serie, seq_length))
      }
      
      begin <- indice
      end <- indice
      sampOut <- currentSampNumber
      currentSampNumber <- sampNumber
      
      ### define the 6 series with higher count for this sample
      tempoDf <- cbind(sum, sampleDf$serie)
      serieSum <- tapply(X = tempoDf[, 1], INDEX = tempoDf[, 2], FUN = sum) 
      serieListe <- sort(unique(sampleDf$serie))
      sumDf <- cbind(serieSum, serieListe)
      sumDf <- sumDf[order(sumDf[, 1], decreasing = TRUE), ]
      sumDf <- as.data.frame(sumDf)
      seqSeries <- sumDf$serieListe[1:6]
      
      # extract data from sampleDf for these 6 series
      table1 <- sampleDf[is.element(sampleDf$serie, seqSeries), ] 
      
      ### loop over the 6 series
      
      # def of sampRepet values
      sampleDfColNames <- colnames(sampleDf)
      
      sampRepetNumbers <- c()
      
      for (n in 1:(length(sampleDfColNames) - 2)) {
        sampRepetNumbers <- c(sampRepetNumbers, 
                              strtoi(substr(sampleDfColNames[n], 
                                            nchar(sampleDfColNames[n]), 
                                            nchar(sampleDfColNames[n]))
                              )
        )
      }
      
      # Iterate over all columns, which holds data for read counts with repeats in columns.
      for (k in 1:(length(sampRepetNumbers))) {
        AlleleDf <- NULL
        genoAllele <- c()
        genoCount <- c()
        scount <- table1[, k] # P1, P2, P3...
        # bind together count for P1, P2... and serie (sampRepet + 1) and seq_length (sampRepet + 2)
        currentTable <- cbind(scount, table1[, (sampRepet + 1):(sampRepet + 2)])
        
        # Do this for top 6 series by count.
        for (seqType in seqSeries) {
          # prepare subtable and store row numbers to complete table1 with results
          rowNbSeqType <- which(table1$serie == seqType)
          if (length(rowNbSeqType != 1)) { # BUG, this always evaluates to integer, which is TRUE -- should be length(x) != 1
            subTable <- currentTable[currentTable$serie == seqType, ]
            ok <- subTable$seq_length # length of sequences for a given series
            
            sv <- subTable[, 1] # select column scount
            out <- rep(0, length(sv))
            i <- 1
            while (i < length(ok)) { # iterate over all but last value
              v0 <- sv[i]
              v0p <- i
              
              if (v0 != 0) {
                ### look for next 3 values if available
                v1 <- 0
                v1p <- 0
                v2 <- 0
                v2p <- 0
                v3 <- 0
                v3p <- 0
                v4 <- 0
                v4p <- 0
                cpt <- i + 1
                
                while (cpt < length(ok) + 1) { # could be cpt <= length(ok)
                  if (ok[cpt] - ok[v0p] == motifLength) {
                    v1 <- sv[cpt]
                    v1p <- cpt
                  }
                  if (ok[cpt]-ok[v0p] == motifLength * 2) {
                    v2 <- sv[cpt]
                    v2p <- cpt
                  }
                  if (ok[cpt] - ok[v0p] == motifLength * 3) {
                    v3 <- sv[cpt]
                    v3p <- cpt
                  }
                  if (ok[cpt] - ok[v0p] == motifLength * 4) {
                    v4 <- sv[cpt]
                    v4p <- cpt
                  }
                  cpt <- cpt + 1
                }
                
                if (v1 * 0.25 > v0) {# modif with threshold
                  if (v2 >= v1) {
                    if (v2 * 0.25 > v1)	{ # v1 potential stutter of v2
                      i <- i + 1									
                    } else {
                      # validate v1 and v2 if v3 < v2 
                      if (v3 < v2) {
                        peakVal <- ok[v1p] #
                        out[v1p] <- peakVal #
                        peakVal <- ok[v2p] #
                        out[v2p] <- peakVal #
                        i <- v2p + 1 #
                      } else {
                        i <- i + 1
                      }
                    }
                  } else {
                    if (v2 >= 0.15 * v1) {
                      if (v3 > v2) {
                        # validate v1
                        peakVal <- ok[v1p]
                        out[v1p] <- peakVal
                        i <- v1p
                      }	
                      else {
                        #valider v1 et v2
                        if (v3p == 0) {
                          i <- v0p + 3
                        } else {
                          i <- v3p
                        }
                        peakVal <- ok[v1p]
                        out[v1p] <- peakVal
                        peakVal <- ok[v2p]
                        out[v2p] <- peakVal
                      }
                    } else {
                      #valider v1
                      i <- v1p
                      ## essai suppresion seul count
                      peakVal <- ok[v1p]
                      out[v1p] <- peakVal
                    }
                  }		
                } else {
                  i <- i + 1
                }
              } else {
                i <- i + 1
              }
            }
            # store data about this serie
            posgeno <- which(out != 0)
            if (length(posgeno != 0)) {
              geno <- paste(out[posgeno], '_', seqType, sep="")
              genoAllele <- c(genoAllele, geno)
              genoCount <- c(genoCount, sv[posgeno])						
            }
          }
        }
        
        # AlleleDF copy in table2
        genoDf <- NULL
        if (length(genoAllele) != 0) {
          genoDf <- data.frame(genoAllele, genoCount)
          genoDf <- genoDf[order(genoDf[, 2], decreasing = TRUE), ]
        }
        if (length(genoAllele) !=0 ) {
          if (dim(genoDf)[1] <= 6) {
            if (dim(genoDf)[1] != 0)
              for (j in 1:(dim(genoDf)[1])) {
                # this does not work if more samples are pooled into one (such as POS)
                targetRowNb <- ((sampOut - 1) * (nrep + 1)) + sampRepetNumbers[k]
                if (nrow(table2) > 846) browser()
                table2[targetRowNb,j] <- as.character(genoDf[j, 1])
                table2[targetRowNb, j + 6] <- genoDf[j, 2]				
              }
          } else {
            for (j in 1:6) {
              targetRowNb <- ((sampOut - 1) * (nrep + 1)) + sampRepetNumbers[k]
              table2[targetRowNb,j] <- as.character(genoDf[j, 1])
              table2[targetRowNb,j + 6] <- genoDf[j, 2]
            }
          }
        }
        
      }
    } else {
      end <- indice
    }
    indice <- indice + 1
  }
  browser()
  fileoutname2 <- paste("table2n_", filename, sep = "")
  ### write table2 to csv file
  write.table(table2, fileoutname2, quote = FALSE, sep = "\t")
}
