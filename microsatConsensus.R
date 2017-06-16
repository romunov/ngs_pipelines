microsatConsensus <- function(filename, repeats, NbEch, threshold, homozygousThreshold) {
  #' Define and add consensus genotype according to the available replicates.
  #' 
  #' @param filename Passed in filename, can be direct or relative path.
  #' @param repeats Number of repeats of each sample.
  #' @param NbEch Number of samples.
  #' @param threshold A threshold to be applied  to sequence count to filter out low quality reads.
  #' @param homozygousThreshold The homozygous threshold.
  #' 
  #' @note All files should be in one folder. To keep things tidy, consider createing a \code{./data} folder or
  #' something similar.
  
  ### script to define and add consensus genotype according to the available replicates
  
  repet <- repeats
  
  ### read  data file 
  data <- read.table(filename, head = TRUE, stringsAsFactors = FALSE)
  
  for (i in 1:NbEch) {
    ### get allele list for each replicate, filtered according to the fixed threshold
    AlleleList=list(c())
    for (j in 1:repet) {
      if (is.na(data[((i-1)*(repet+1))+j,1])) {
        AlleleList[j]=list(c(NA))
      } else {
        liste=list(c())
        #liste=list()
        listCount=list(c())
        #listCount=list()
        liste[1]=data[((i-1)*(repet+1))+j,1]
        listCount[1]=data[((i-1)*(repet+1))+j,7]
        cpt=8
        go=1
        while (cpt<13 & go==1)
        {
          if (!is.na(data[((i-1)*(repet+1))+j,cpt]))
          {
            if (data[((i-1)*(repet+1))+j,cpt]>=(threshold*data[((i-1)*(repet+1))+j,7]))
            {
              # get the size of the current allele
              splitname=strsplit(data[((i-1)*(repet+1))+j,cpt-6],'_')
              currentAlleleLength=as.numeric(splitname[[1]][1])						
              
              
              go2=1
              # check for allele of the same size
              listLength=length(liste)
              for (l in 1:listLength)
              {
                splitname=strsplit(liste[[l]][1],'_')
                AlleleLength=as.numeric(splitname[[1]][1])
                
                if (currentAlleleLength==AlleleLength)
                {
                  go2=0
                  if (data[((i-1)*(repet+1))+j,cpt]>0.3*listCount[[l]][1])
                  {
                    liste[listLength+1]=data[((i-1)*(repet+1))+j,cpt-6]
                    listCount[listLength+1]=data[((i-1)*(repet+1))+j,cpt]
                  }
                }
                else
                {
                  #liste[listLength+1]=data[((i-1)*(repet+1))+j,cpt-6]
                  #listCount[listLength+1]=data[((i-1)*(repet+1))+j,cpt]
                }
                
                
                
              }
              
              if (go2==1)
              {
                liste[listLength+1]=data[((i-1)*(repet+1))+j,cpt-6]
                listCount[listLength+1]=data[((i-1)*(repet+1))+j,cpt]
              }
              
              cpt=cpt+1
              
            }
            else
            {
              go=0
            }
          }
          else
          {
            go=0
          }
        }
        
        AlleleList[[j]]=liste
      }
      
    }	
    
    ### concat allele list
    Alleles=unlist(AlleleList)
    tabAlleles=sort(table(Alleles),decreasing=TRUE)
    
    if (length(tabAlleles)==0)
    {
      data[((i-1)*(repet+1))+(repet+1),1]=0
      data[((i-1)*(repet+1))+(repet+1),2]=0
    } 
    else
    {
      ### select only alleles with count>=2
      tabAllelesSelect=tabAlleles[tabAlleles>=2]
      
      if (length(tabAllelesSelect)==0)
      {
        data[((i-1)*(repet+1))+(repet+1),1]=0
        data[((i-1)*(repet+1))+(repet+1),2]=0
      }
      
      if (length(tabAllelesSelect)==1) 
      {
        if (tabAllelesSelect[1]>=homozygousThreshold*repet) 
        {
          data[((i-1)*(repet+1))+(repet+1),1]=names(tabAllelesSelect)
          data[((i-1)*(repet+1))+(repet+1),2]=names(tabAllelesSelect)
          # quality
          data[((i-1)*(repet+1))+(repet+1),3]=tabAllelesSelect[1]/repet
          
          
        }
        else
        {
          data[((i-1)*(repet+1))+(repet+1),1]=0
          data[((i-1)*(repet+1))+(repet+1),2]=0
        }
      }
      
      if (length(tabAllelesSelect)==2)
      {
        all1=names(tabAllelesSelect)[1]
        all2=names(tabAllelesSelect)[2]
        len_all1=as.numeric(strsplit(all1,'_')[[1]][1])
        len_all2=as.numeric(strsplit(all2,'_')[[1]][1])
        if (len_all2<len_all1)
        {
          data[((i-1)*(repet+1))+(repet+1),1]=all2
          data[((i-1)*(repet+1))+(repet+1),2]=all1
        }
        else 
        {
          data[((i-1)*(repet+1))+(repet+1),1]=all1
          data[((i-1)*(repet+1))+(repet+1),2]=all2
        }
        
        
        #data[((i-1)*(repet+1))+(repet+1),1]=names(tabAllelesSelect)[1]
        #data[((i-1)*(repet+1))+(repet+1),2]=names(tabAllelesSelect)[2]
        
        # quality
        allelesNames=names(tabAllelesSelect)
        genoCount=0
        for (n in 1:repet)
        {
          if (is.element(allelesNames[1],AlleleList[[n]]) & is.element(allelesNames[2],AlleleList[[n]]))
          {
            genoCount=genoCount+1
          }		
        }
        data[((i-1)*(repet+1))+(repet+1),3]=genoCount/repet
      }
      
      if (length(tabAllelesSelect)>2) {
        data[((i-1)*(repet+1))+(repet+1),1]=1
        data[((i-1)*(repet+1))+(repet+1),2]=1
      }
      
    }
  }
  
  ### print updated table
  fileoutname <- paste("c", filename, sep = "")
  write.table(data, filename, quote = FALSE)
  
  ### create table_cns
  table_cns <- data[seq(repet + 1, NbEch * (repet + 1), (repet + 1)), ]
  table_cns <- table_cns[, 1:3]
  colnames(table_cns) <- c("allele1", "allele2", "quality")
  fileoutname2 <- paste("cns_", filename, sep = "")
  write.table(table_cns, fileoutname2, quote = FALSE)
  
  ### check for ref file
  # parts=strsplit(filename,'_')
  # [[1]]
  # [1] "table2n"      "MICROSAT.PCR" "M1P5P.MM"     "25"          
  # [5] "serie.tab"   
  # reffilename=paste(substr(parts[[1]][3],1,2),'_',parts[[1]][4],'_ref.txt',sep="")
  # [1] "M1_25_ref.txt"
  
  # extract bits and pieces from the filename
  samplename <- gsub("(^table2n_MICROSAT\\.PCR_)(.*_\\d+)(_serie\\.tab$)", "\\2", filename) # includes locus name
  locusname <- gsub("(^table2n_MICROSAT\\.PCR_)(.*_)(\\d+)(_serie\\.tab$)", "\\3", filename)
  sn <- gsub("(^table2n_MICROSAT\\.PCR_)(.*)_(\\d+)_(serie\\.tab$)", "\\2", filename)
  reffilename <- paste(paste("UA_", locusname, sep = ""), "_ref.txt", sep = "")
  
  # find serie file
  find.file <- paste("^MICROSAT.PCR_", samplename, "_serie.tab$", sep = "")
  locate.file <- list.files(pattern = find.file)
  
  if (file.exists(reffilename)) {
    #### modify cns table allele names according to ref file
    # read ref file
    # parts=strsplit(filename,'_')
    # reffilename=paste(substr(parts[[1]][3],1,2),'_',parts[[1]][4],'_ref.txt',sep="")
    
    refAlleles <- read.table(reffilename, header = TRUE, stringsAsFactors = FALSE)
    
    # create filename to store allele name changes
    # allele_changes_filename <- paste(parts[[1]][3],'_',parts[[1]][4],"_allele_changes.txt",sep="")
    allele_changes_filename <- paste(samplename, "_allele_change.txt", sep = "")
    # create new file every time locus is being processed
    cat("sample\toldname\tnewname", file = allele_changes_filename, sep = "\n")
    
    # read data serie file
    dataserie <- read.table(locate.file, header = T)
    
    # parse cns table and check each allele
    table_cnsh <- table_cns
    # reads through the genotype table and tries to update allele names on the fly
    for (i in 1:NbEch) {
      for (j in 1:2) {
        allele <- table_cnsh[i, j]
        if (allele != '0' & allele != '1') {
          size <- as.numeric(strsplit(allele, '_')[[1]][1])
          serie <- as.numeric(strsplit(allele, '_')[[1]][2])
          sub <- dataserie[dataserie$seq_length == size & dataserie$series == serie, ]
          seq <- as.character(sub$sequence)
          found <- 0
          namefound <- 0
          cpt <- 1
          while (cpt <= (dim(refAlleles)[1]) & found == 0) {
            if (seq == refAlleles[cpt, 2]) {
              found <- 1
              if (allele != refAlleles[cpt, 1]) { # if new sequence found
                table_cnsh[i,j] <- refAlleles[cpt, 1]
                browser()
                # ... write to allele_change_filename
                cat(rownames(table_cnsh)[i], allele, refAlleles[cpt, 1], "\n", 
                    file = allele_changes_filename, sep = "\t", append = TRUE)
              }
            }
            cpt <- cpt + 1
          }
          if (found == 0) {
            refAllelesNames <- refAlleles[, 1]
            if (!is.element(allele, refAllelesNames)) {
              refAlleles <- rbind(refAlleles, c(allele, seq, samplename = sn, locus = locusname))
            } else {
              ct <- 1
              new_name <- paste(size, '_', strtoi(ct), sep = '')
              while (is.element(new_name, refAllelesNames)) {
                ct = ct + 1
                new_name = paste(size, '_', strtoi(ct), sep = '')
              }
              refAlleles <- rbind(refAlleles, c(new_name, seq, samplename = sn, locus = locusname))
              table_cnsh[i, j] <- new_name			
            }
          }
        }
      }
    }
    ### write table_cnsh
    fileoutnamecnsh <- paste("cnsh_", filename, sep = "")
    write.table(table_cnsh, fileoutnamecnsh, quote = FALSE)
    
    ### write updated table refAlleles
    write.table(refAlleles, reffilename, quote = FALSE, row.names = FALSE)
    
  }
  if (!file.exists(reffilename)) {
    ### create allele ref file
    refAlleles <- c(table_cns$allele1, table_cns$allele2)
    refAlleles <- unique(refAlleles)
    refAlleles <- refAlleles[refAlleles != 1 & refAlleles != 0]
    
    if (length(locate.file) == 0) {
      stop(sprintf("Unable to locate %s, please check that the file is in the same folder as input.", find.file))
    }
    
    dataserie <- read.table(locate.file, header = TRUE)
    
    if (length(refAlleles) > 0) {
      # This chunk has been disabled because it's having problems when no alleles are present.
      seqs <- c()
      
      for (i in 1:(length(refAlleles))) {
        size <- as.numeric(strsplit(refAlleles[i], '_')[[1]][1])
        serie <- as.numeric(strsplit(refAlleles[i], '_')[[1]][2])
        sub <- dataserie[dataserie$seq_length == size & dataserie$series == serie, ]
        seqs <- c(seqs, as.character(sub$sequence))
      }
      
      dfRefAlleles <- data.frame(refAlleles, seqs, samplename = sn, locus = locusname)
      
      ### write dfRefAlleles to ref file
      # tmp=strsplit(filename,'_')
      # reffilename=paste(substr(tmp[[1]][3],1,2),'_',strtoi(tmp[[1]][4]),"_ref.txt",sep="")
      write.table(dfRefAlleles, reffilename, quote = FALSE, row.names = FALSE)
      
      ### write table_cnsh
      fileoutnamecnsh=paste("cnsh_", filename, sep = "")
      write.table(table_cns, fileoutnamecnsh, quote = FALSE)
    } else {
      message(sprintf("Skipping %s because no alleles found.", paste(sn, "_", locusname, sep = "")))
    }
  }
}
