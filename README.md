This is how I analyzed GATC data for the `DinAlp bear` project.

```
# extract files
gzip -d NG-12425_E45175_lib187741_5356_6_1.fastq.gz
gzip -d NG-12425_E45175_lib187741_5356_6_2.fastq.gz

# move them to wherever you will be running the pipeline 
mv * ~/Documents/gatc
# also add .ngsfilter to same position where .fastq files are

# combine and make consensus reads from tw 1D reads
illuminapairedend -r NG-12425_E45175_lib187741_5356_6_2.fastq NG-12425_E45175_lib187741_5356_6_1.fastq | tee rawdata_gatc_medvedi.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_gatc_medvedi.
 
# close files
touch rawdata_gatc_medvedi.Bad.fastq
touch rawdata_gatc_medvedi.Alignement.fastq
touch rawdata_gatc_medvedi.fastq
 
# recognize and assign locus and sample name
ngsfilter -t UA_gatc_8_plate.ngsfilter -u rawdata_gatc_medvedi.unidentified.fastq  rawdata_gatc_medvedi.Alignement.fastq > rawdata_gatc_medvedi.filtered.fastq
 
# cut by locus
obisplit -p MICROSAT.PCR_ -t experiment rawdata_gatc_medvedi.filtered.fastq
 
# find unique sequences and count them
obiuniq -m sample MS.PCR_UA_MxRout1_03.fastq > MS.PCR_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_06.fastq > MS.PCR_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_14.fastq > MS.PCR_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_16.fastq > MS.PCR_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_17.fastq > MS.PCR_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_25.fastq > MS.PCR_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_51.fastq > MS.PCR_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_57.fastq > MS.PCR_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_63.fastq > MS.PCR_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_64.fastq > MS.PCR_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_65.fastq > MS.PCR_UA_MxRout1_65.uniq.fasta & 
obiuniq -m sample MS.PCR_UA_MxRout1_67.fastq > MS.PCR_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_68.fastq > MS.PCR_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MS.PCR_UA_MxRout1_ZF.fastq > MS.PCR_UA_MxRout1_ZF.uniq.fasta
 
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_03.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_06.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_14.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_16.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_17.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_25.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_51.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_57.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_63.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_64.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_65.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_67.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_68.uniq.tab
obigrep -p 'count>1' MICROSAT.PCR_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_UA_MxRout1_ZF.uniq.tab

# first edit the script for the number of used cores. Result is genotypes_UA_GATC.txt file.
Rscript R_processing.R
```
