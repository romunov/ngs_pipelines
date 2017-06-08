### script to analyse microsat NGS tab file per sample through sequence family
microsatAlleles <- function(nrep, motifLength, numsamples, numseqfamily, input) {
  nrep <- as.numeric(nrep)
  motifLength <- as.numeric(motifLength)
  NbEch <- as.numeric(numsamples)
  minNbOfSeqPerFamily <- as.numeric(numseqfamily)
  filename <- input
  
  ### read csv data file and remove def column
  #data=read.csv(filename,h=T,sep="\t")
  data <- read.table(filename, head = TRUE)
  
  ######### filtering according to the minimal number of sequences per family provided
  familySeqCount=tapply(data$count, data$series, sum)
  series=unique(data$series)
  familydf=as.data.frame(cbind(familySeqCount,series))
  seriesok=familydf$series[familydf$familySeqCount>minNbOfSeqPerFamily]
  data=data[is.element(data$series,seriesok),]
  
  
  
  ### def first and last sample positions (sample...R.) (voir regexp)
  # def repet tags to use according to repet number
  Rrep <- paste("P", 1:nrep, sep = "")
  
  # sort data according to seq length values
  data=data[order(data$seq_length),]
  
  ### get sample numbers and experiment name
  # df=data[,posFirstSample:posLastSample] 
  df <- data[, colnames(data)[grepl("_P\\d+$", colnames(data))]]
  
  nom <- colnames(df)
  dataColNb <- dim(data)[2]
  ext <- as.numeric(sapply(nom, function(x) strsplit(x,'_')[[1]][2]))
  sampleNb <- unique(ext)
  nExt <- strsplit(nom[1],'_')[[1]][1]
  nameExp <- substr(nExt, 8, nchar(nExt))
  
  
  ###### prepare table2 to transfer genotype data id + 6 geno + 6 counts
  ### create and initialize table 2
  table2=NULL
  table2=as.data.frame(matrix(data=NA,nrow=NbEch*(nrep+1),ncol=12))
  
  namesTable2=c()
  for (i in 1:(NbEch)) {
    samp=i
    add=""
    if (samp<10) 
    {
      add="00"
    }
    if (samp>9 & samp <100)
    {
      add="0"
    }
    
    for (j in 1:(nrep))
    {
      currentName=paste("sample.",nameExp,'_',add,strtoi(samp),'_R',strtoi(j),sep="")
      namesTable2=c(namesTable2,currentName)
    }
    cnsName=paste("sample.",nameExp,'_',add,strtoi(samp),'_CNS',sep="")
    namesTable2=c(namesTable2,cnsName)
  }
  
  rownames(table2) <- namesTable2
  colNamesTable2 <- c("allele1","allele2","allele3","allele4","allele5","allele6","count_allele1","count_allele2","count_allele3","count_allele4","count_allele5","count_allele6")
  colnames(table2) <- colNamesTable2
  
  ### loop over sample columns
  dfColNames=colnames(df)
  indice=1
  begin=1
  end=1
  # init sampNumber
  currentSampNumber=as.numeric(strsplit(dfColNames[1],'_')[[1]][2])
  
  while (indice<=(length(dfColNames)+1)) {
    if (indice<=(length(dfColNames)))
    {
      sampNumber=as.numeric(strsplit(dfColNames[indice],'_')[[1]][2])
    }
    else 
    {
      sampNumber="" ###  to analyze the last sample
    }
    if (currentSampNumber!=sampNumber & indice!=1) ### ajouter ici une condition pour traiter le dernier echantillon ou==length()
    {
      
      ### extract corresponding data
      sampleDf=df[,(begin):(end)]
      if (begin==end)
      {
        
        sampRepet=1
        sum=sampleDf
        serie=data$serie
        seq_length=data$seq_length
        sample=sampleDf
        sampleDf=as.data.frame(cbind(sample,serie,seq_length))
        namesCol=c(dfColNames[indice-1],"serie","seq_length")
        colnames(sampleDf)=namesCol
      }
      else
      {
        sampRepet=dim(sampleDf)[2]
        sum=apply(sampleDf,1,sum)
        serie=data$serie
        seq_length=data$seq_length
        sampleDf=as.data.frame(cbind(sampleDf,serie,seq_length))
      }
      
      begin=indice
      end=indice
      sampOut=currentSampNumber
      currentSampNumber=sampNumber
      
      
      ### define the 6 series with higher count for this sample
      tempoDf=cbind(sum,sampleDf$serie)
      serieSum=tapply(tempoDf[,1],tempoDf[,2],sum) 
      serieListe=sort(unique(sampleDf$serie))
      sumDf=cbind(serieSum,serieListe)
      sumDf=sumDf[order(sumDf[,1],decreasing=TRUE),]
      sumDf=as.data.frame(sumDf)
      seqSeries=sumDf$serieListe[1:6]
      
      # extract data from sampleDf for these 6 series
      table1=sampleDf[is.element(sampleDf$serie,seqSeries),] 
      
      ### loop over the 6 series
      
      # def of sampRepet values
      sampleDfColNames=colnames(sampleDf)
      
      sampRepetNumbers=c()
      for (n in 1:(length(sampleDfColNames)-2))
      {
        
        sampRepetNumbers=c(sampRepetNumbers,strtoi(substr(sampleDfColNames[n],nchar(sampleDfColNames[n]),nchar(sampleDfColNames[n]))))
      }
      
      for (k in 1:(length(sampRepetNumbers)))
        
      {	
        AlleleDf=NULL
        genoAllele=c()
        genoCount=c()
        scount=table1[,k]
        currentTable=cbind(scount,table1[,(sampRepet+1):(sampRepet+2)])
        
        for (seqType in seqSeries) 
        {
          # prepare subtable and store row numbers to complete table1 with results
          rowNbSeqType=which(table1$serie==seqType)
          if (length(rowNbSeqType!=1))
          {
            subTable=currentTable[currentTable$serie==seqType,]
            ok=subTable$seq_length
            
            sv=subTable[,1]
            out=rep(0,length(sv))
            i=1
            while (i<length(ok))
            {
              v0=sv[i]
              v0p=i
              
              if (v0!=0)
                ### look for next 3 values if available
              {
                v1=0
                v1p=0
                v2=0
                v2p=0
                v3=0
                v3p=0
                v4=0
                v4p=0
                cpt=i+1
                
                while (cpt<length(ok)+1)
                {
                  if (ok[cpt]-ok[v0p]==motifLength)
                  {
                    v1=sv[cpt]
                    v1p=cpt
                  }
                  if (ok[cpt]-ok[v0p]==motifLength*2)
                  {
                    v2=sv[cpt]
                    v2p=cpt
                  }
                  if (ok[cpt]-ok[v0p]==motifLength*3)
                  {
                    v3=sv[cpt]
                    v3p=cpt
                  }
                  if (ok[cpt]-ok[v0p]==motifLength*4)
                  {
                    v4=sv[cpt]
                    v4p=cpt
                  }
                  cpt=cpt+1
                }
                
                if (v1*0.25>v0) # modif with threshold
                {
                  if (v2>=v1) # modif >=
                  {
                    
                    if (v2*0.25>v1)	#v1 potential stutter of v2
                    {
                      i=i+1		#									
                    }
                    else #
                    {
                      # validate v1 and v2 if v3<v2 
                      if (v3<v2) #
                      {
                        peakVal=ok[v1p] #
                        out[v1p]=peakVal #
                        peakVal=ok[v2p] #
                        out[v2p]=peakVal #
                        i=v2p+1 #
                      }
                      else
                      {
                        i=i+1
                      }
                    }
                    
                    
                  }
                  else
                  {
                    if (v2>=0.15*v1)
                    {
                      if (v3>v2) 
                      {
                        
                        # validate v1
                        peakVal=ok[v1p]
                        out[v1p]=peakVal
                        i=v1p
                        
                        
                      }	
                      else
                      {
                        #valider v1 et v2
                        
                        if (v3p==0) 
                        {
                          i=v0p+3
                        }
                        else
                        {
                          i=v3p
                        }
                        
                        
                        peakVal=ok[v1p]
                        
                        out[v1p]=peakVal
                        
                        
                        peakVal=ok[v2p]
                        
                        out[v2p]=peakVal
                        
                      }
                      
                    }
                    else
                    {
                      #valider v1
                      i=v1p
                      
                      ## essai suppresion seul count
                      peakVal=ok[v1p]
                      ##
                      
                      
                      out[v1p]=peakVal
                    }
                  }		
                }
                else
                {
                  i=i+1
                }
              }
              else
              {
                i=i+1
              }
            }
            
            
            # store data about this serie
            posgeno=which(out!=0)
            if (length(posgeno!=0))
            {
              geno=paste(out[posgeno],'_',seqType,sep="")
              genoAllele=c(genoAllele,geno)
              genoCount=c(genoCount,sv[posgeno])						
            }
            
            
            
          }
        }
        
        # AlleleDF copy in table2
        genoDf=NULL
        if (length(genoAllele)!=0)
        {
          genoDf=data.frame(genoAllele,genoCount)
          genoDf=genoDf[order(genoDf[,2],decreasing=TRUE),]
        }
        if (length(genoAllele)!=0)
        {
          if (dim(genoDf)[1]<=6)
          {
            if (dim(genoDf)[1]!=0)
              for (j in 1:(dim(genoDf)[1]))
              {
                targetRowNb=((sampOut-1)*(nrep+1))+sampRepetNumbers[k]
                table2[targetRowNb,j]=as.character(genoDf[j,1])
                table2[targetRowNb,j+6]=genoDf[j,2]				
              }
            
          }
          else
          {
            for (j in 1:6)
            {
              targetRowNb=((sampOut-1)*(nrep+1))+sampRepetNumbers[k]
              table2[targetRowNb,j]=as.character(genoDf[j,1])
              table2[targetRowNb,j+6]=genoDf[j,2]
            }
            
          }
        }
        
      }
      
    }
    
    else
      
    {
      end=indice
    }
    indice=indice+1
    
  }
  
  fileoutname2=paste("table2n_",filename,sep="")
  ### write table2 to csv file
  write.table(table2,fileoutname2,quote=FALSE,sep="\t")
  
}

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

microsatConsensus <- function(filename, repeats, NbEch, threshold, homozygousThreshold) {
  ### script to define and add consensus genotype according to the available replicates
  
  #### collect and check arguments provided
  # args <- commandArgs(TRUE)
  # 
  # if (length(args)<5 | length(args)>5)
  # {
  #   print ("##### Error Message #####")
  #   print ("you need to provide five arguments in this order :")
  #   print ("1- The number of repeats for each sample")
  #   print ("2- The number of samples")
  #   print ("3- The threshold to apply to seq count")
  #   print ("4- The homozygous threshold")
  #   print ("5- The name of the input tab file")	
  #   stop()
  # }
  
  repet <- repeats
  
  ### read  data file 
  data=read.table(filename,h=T,stringsAsFactors=FALSE)
  
  for (i in 1:NbEch) {
    ### get allele list for each replicate, filtered according to the fixed threshold
    AlleleList=list(c())
    for (j in 1:repet) 
    {
      #print (j)
      if (is.na(data[((i-1)*(repet+1))+j,1])) 
      {
        AlleleList[j]=list(c(NA))
      }
      else
      {
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
      
      if (length(tabAllelesSelect)>2)
      {
        data[((i-1)*(repet+1))+(repet+1),1]=1
        data[((i-1)*(repet+1))+(repet+1),2]=1
      }
      
    }
  }
  
  
  ### print updated table
  fileoutname=paste("c",filename,sep="")
  write.table(data,fileoutname,quote=FALSE)
  
  ### create table_cns
  table_cns=data[seq(repet+1,NbEch*(repet+1),(repet+1)),]
  table_cns=table_cns[,1:3]
  colnames(table_cns)=c("allele1","allele2","quality")
  fileoutname2=paste("cns_",filename,sep="")
  write.table(table_cns,fileoutname2,quote=FALSE)
  
  
  ### check for ref file
  samplename <- gsub("(^table2n_MICROSAT\\.PCR_)(.*_\\d+)(_serie\\.tab$)", "\\2", filename)
  reffilename <- paste(samplename, "_ref.txt", sep = "")
  
  if (file.exists(reffilename)) {
    #### modify cns table allele names according to ref file
    # read ref file
    #parts=strsplit(filename,'_')
    #reffilename=paste(substr(parts[[1]][3],1,2),'_',parts[[1]][4],'_ref.txt',sep="")
    refAlleles=read.table(reffilename,h=T,stringsAsFactors=FALSE)
    
    # create filename to store allele name changes
    allele_changes_filename=paste(parts[[1]][3],'_',parts[[1]][4],"_allele_changes.txt",sep="")
    cat("sample\toldname\tnewname",file=allele_changes_filename,sep="\n")
    
    # read data serie file
    parts=strsplit(filename,'_')
    seriefilename=paste(parts[[1]][2],'_',parts[[1]][3],'_',parts[[1]][4],'_',parts[[1]][5],sep="")
    dataserie=read.table(seriefilename,h=T)
    
    
    # parse cns table and check each allele
    table_cnsh=table_cns
    for (i in 1:NbEch) {
      for (j in 1:2) {
        #print (i)
        #print (j)
        allele=table_cnsh[i,j]
        #print (allele)
        if (allele!='0' & allele!='1')
        {
          size=as.numeric(strsplit(allele,'_')[[1]][1])
          serie=as.numeric(strsplit(allele,'_')[[1]][2])
          sub=dataserie[dataserie$seq_length==size & dataserie$series==serie,]
          seq=as.character(sub$sequence)
          found=0
          namefound=0
          cpt=1
          while (cpt<=(dim(refAlleles)[1]) & found==0)
          {
            if (seq==refAlleles[cpt,2])
            {
              found=1
              if (allele!=refAlleles[cpt,1])
              {
                table_cnsh[i,j]=refAlleles[cpt,1]
                cat(rownames(table_cnsh)[i],allele,refAlleles[cpt,1],"\n",file=allele_changes_filename,sep="\t",append=TRUE)
              }
            }
            cpt=cpt+1
          }
          if (found==0)
          {
            refAllelesNames=refAlleles[,1]
            if (!is.element(allele,refAllelesNames))
            {
              refAlleles=rbind(refAlleles,c(allele,seq))
            }
            else
            {
              ct=1
              new_name=paste(size,'_',strtoi(ct),sep='')
              while (is.element(new_name,refAllelesNames))
              {
                ct=ct+1
                new_name=paste(size,'_',strtoi(ct),sep='')
              }
              refAlleles=rbind(refAlleles,c(new_name,seq))
              table_cnsh[i,j]=new_name			
            }
            
          }
        }
      }
      
    }
    ### write table_cnsh
    fileoutnamecnsh=paste("cnsh_",filename,sep="")
    write.table(table_cnsh,fileoutnamecnsh,quote=FALSE)
    
    ### write updated table refAlleles
    write.table(refAlleles,reffilename,quote=FALSE,row.names=FALSE)
    
  }
  if (!file.exists(reffilename)) {
    ### create allele ref file
    refAlleles=c(table_cns$allele1,table_cns$allele2)
    refAlleles=unique(refAlleles)
    refAlleles=refAlleles[refAlleles!=1 & refAlleles!=0]
    
    ### read serie file
    # parts=strsplit(filename,'_')
    # seriefilename=paste(parts[[1]][2],'_',parts[[1]][3],'_',parts[[1]][4],'_',parts[[1]][5],sep="")
    
    find.file <- paste("^MICROSAT.PCR_", samplename, "_serie.tab$", sep = "")
    locate.file <- list.files(pattern = find.file)
    
    if (length(locate.file) == 0) {
      stop(sprintf("Unable to locate %s, please check that the file is in the same folder as input.", find.file))
    }
    
    dataserie <- read.table(locate.file, header = TRUE)
    
    if (FALSE) {
      # This chunk has been disabled because it's having problems when no alleles are present.
      seqs=c()
      
      for (i in 1:(length(refAlleles))) {
        size=as.numeric(strsplit(refAlleles[i],'_')[[1]][1])
        serie=as.numeric(strsplit(refAlleles[i],'_')[[1]][2])
        sub=dataserie[dataserie$seq_length==size & dataserie$series==serie,]
        seqs=c(seqs,as.character(sub$sequence))
      }
      
      dfRefAlleles=data.frame(refAlleles,seqs)
      
      ### write dfRefAlleles to ref file
      # samplename <- gsub("(^table2n_MICROSAT\\.PCR_)(.*_\\d+)(_serie\\.tab$)", "\\2", filename)
      # reffilename <- paste(samplename, "_ref.txt", sep = "")
      write.table(dfRefAlleles,reffilename,quote=FALSE,row.names=FALSE)
      
      ### write table_cnsh
      fileoutnamecnsh=paste("cnsh_",filename,sep="")
      write.table(table_cns,fileoutnamecnsh,quote=FALSE)
    }
  }
  
  
  
  
}

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
  
  # browser()
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
