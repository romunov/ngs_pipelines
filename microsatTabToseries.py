#!/usr/bin/python

"""  add seq family data as an additional column to microsat tab file """

import sys
from optparse import OptionParser
import re

### argument parsing
parser=OptionParser("Usage : %prog arg1 arg2")
parser.add_option("-f", "--filename", dest="filename", type="string", help="Name of the tab file")
parser.add_option("-m", "--motif", dest="motif", type="string", help="DNA sequence of motif")
(options, args) = parser.parse_args()

if len(sys.argv[1:]) !=4 :
    sys.stdout.write("vous devez entrer 2 arguments")
    print ""
    sys.stdout.write("-h pour de l'aide")
    print ""
    sys.exit(0)

def mean (list) :
	cpt=0
	total=0
	for val in list :
		total=total+val
		cpt=cpt+1
	meanVal=total/cpt 
	return meanVal

filename=options.filename
motif=options.motif

lenM=len(motif)
file=open(filename,'r')

## define fileoutname
# parts=filename.split('.')
# User regular expression to find anything before .uniq.tab.
parts = re.search(pattern="(^MICROSAT\\.PCR_.*)(\\.uniq\\.tab)", string=filename)
# fileoutname=parts[0]+'.'+parts[1]+"_serie.tab"
fileoutname = parts.group(1) + "_serie.tab"
fileout = open(fileoutname,'w')

first=1
nextType=1
info="skip"

currentList=[]
for i in range(0,10) :
	currentList.append("")

seq1_all=[]
seq2_all=[]
motifm1_all=[]
posmotifm1_all=[]
motifm2_all=[]
posmotifm2_all=[]
motifm3_all=[]
posmotifm3_all=[]
motifm4_all=[]
posmotifm4_all=[]
	
for line in file :
	if (info=="skip") :
	  # add a column name 'series'
		sortie = line[0:len(line)-1] + "\t" + "series"
		fileout.write(sortie)
		fileout.write('\n')
		info="go"
	else :
		element=line.split()
		seqsat=element[len(element)-1]

		### add a test to remove sequences that do not contain two motifs
		found=0
		cpt=-1
		if len(seqsat)>lenM*2 :
			
			while (found==0 and cpt<(len(seqsat)-lenM*2)) :
				cpt=cpt+1
				seq=seqsat[cpt:cpt+lenM*2]
				if (seq==motif*2) :
					
					found=1
		
		if found==1 :
			### search for seq1 and repeat start
			found=0
			cpt=-1
			while (found==0 and cpt<(len(seqsat)-lenM)) :
				cpt=cpt+1
				seq=seqsat[cpt:cpt+lenM]
				if (seq==motif) :
					found=1
			endseq1=cpt-1
			seq1=seqsat[0:(endseq1+1)]	
			currentList[0]=seq1
	
			### search for seq2 and repeat end
			found=0
			cpt=len(seqsat)-lenM*2+1
			while (found==0 and cpt>-1) :
				cpt=cpt-1
				seq=seqsat[cpt:cpt+lenM*2]
				if (seq==motif*2) :
					found=1
			startseq2=cpt+lenM*2
			seq2=seqsat[startseq2:len(seqsat)]
			currentList[9]=seq2
	
			### search for modified motif
			debm=0
			endm=0
			cpt=endseq1+1
			nbmotif=0
			while (cpt<(startseq2-lenM)) :
				currentseq=seqsat[cpt:cpt+lenM]
				if (currentseq==motif) :
					cpt=cpt+lenM
				else :
					debm=cpt
					while (currentseq!=motif) :
						cpt=cpt+1
						currentseq=seqsat[cpt:cpt+lenM]
					endm=cpt-1
					motifm=seqsat[debm:endm+1]
					nbmotif=nbmotif+1
					### actual limit of 4 modified motifs
					if (nbmotif<5) :
						currentList[(nbmotif*2)-1]=motifm
						currentList[nbmotif*2]=debm
	
			### search into family data
			if (first==1) :
				seq1_all.append(currentList[0])	
				motifm1_all.append(currentList[1])
				posmotifm1_all.append(currentList[2])
				motifm2_all.append(currentList[3])
				posmotifm2_all.append(currentList[4])
				motifm3_all.append(currentList[5])
				posmotifm3_all.append(currentList[6])
				motifm4_all.append(currentList[7])
				posmotifm4_all.append(currentList[8])
				seq2_all.append(currentList[9])
				sortie=line[0:len(line)-1]+"\t"+str(nextType)
				fileout.write(sortie)
				fileout.write('\n')
				nextType=nextType+1
				first=2

	
			else :
				cpt=0
				found=0
				while (found==0 and cpt<len(seq1_all)) :
				
					if (seq1_all[cpt]==currentList[0] and motifm1_all[cpt]==currentList[1] and posmotifm1_all[cpt]==currentList[2] and motifm2_all[cpt]==currentList[3] and posmotifm2_all[cpt]==currentList[4] and motifm3_all[cpt]==currentList[5] and posmotifm3_all[cpt]==currentList[6] and motifm4_all[cpt]==currentList[7] and posmotifm4_all[cpt]==currentList[8] and seq2_all[cpt]==currentList[9]) :
						seqType=cpt+1
						found=1
						sortie=line[0:len(line)-1]+"\t"+str(seqType)
						fileout.write(sortie)
						fileout.write('\n')
					cpt=cpt+1
				if (found==0) :
					seq1_all.append(currentList[0])				
					motifm1_all.append(currentList[1])
					posmotifm1_all.append(currentList[2])
					motifm2_all.append(currentList[3])
					posmotifm2_all.append(currentList[4])
					motifm3_all.append(currentList[5])
					posmotifm3_all.append(currentList[6])
					motifm4_all.append(currentList[7])
					posmotifm4_all.append(currentList[8])
					seq2_all.append(currentList[9])
					seqType=nextType
					nextType=nextType+1
					sortie=line[0:len(line)-1]+"\t"+str(seqType)
					fileout.write(sortie)
					fileout.write('\n')
			
			currentList=[]
			for i in range(0,10) :
				currentList.append("")
	
	
file.close()




    	
    


    
    



    
   
  

    
    
