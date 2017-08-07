This is how I analyzed NGS data for the `DinAlp bear` project.

1) Create NGS filters. Help yourself with script `batch_create_ngs_filters.R`. Make sure you designate correct folders for aliquot plates.
2) Run the obitools pipeline below. It helps if you have a sufficient number of cores at your disposal. :)
3) Once you get the `series.tab` files you will be able to process the genotypes in R using `library(fishbone)`* designed for allele calling. See `NGS_pipeline.R`.

* You can install this using `devtools::install_github("romunov/fishbone")`, but you will need the appropriate tool chain.

Files/folders in `/DAB` are:
* `pars.csv` - list of parameters used for allele calling, see function `callAllele` for all the information
* `/1_ngsfilters_hiseq1`: ngs filters for hiseq run 1
* `/2_uniq_tab_hiseq1`: result of obitools pipeline
* `/3_lib_sample_locus`: final result of obitools pipeline, incl. \*serie.tab files produced by `microsatTabToseries.py` script
* `/aliquot_plates`: which plates comprise a library and which sample was put into which position
* `/data`: this is where the genotypes are at before they're shipped off

To determine genetic sex, use script `determine_sex_from_ZF.R` on `*ZF.uniq.tab`. Make sure you correctly set parameters in the INPUT section. Briefly, the algorithm checks for disbalance of X and Y alleles and calles a male or female appropriately. See algorithm details in script comments.

This is the pipeline used to download and sift through the data. Mind you that intermediate files can grow large and may need to be moved or deleted before running the next step.

```
Fasteris DinAlpBear HiSeq2
# aliquot plate od 45-52 sem naredil na roke (Bine in Maja sta bila takrat na dopustu) - možen vir napak

wget -r https://data.fasteris.com/private/JFV/JFV-13-24/data/ --no-parent --user=tomaz.skrbinsek@gmail.com --password=k8zS3sh1 -A fastq.gz --no-check-certificate

# move .gz files to archive in /mnt/ngs_data_storage

# This will remove all extracted .gz files, so make sure you make a copy of them before proceeding
ls *.gz | parallel gunzip

# process this on SSD
illuminapairedend -r 170721_SND405_A_L001_JFV-13_R2.fastq 170721_SND405_A_L001_JFV-13_R1.fastq | tee rawdata_L001_JFV-13_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-13_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-14_R2.fastq 170721_SND405_A_L001_JFV-14_R1.fastq | tee rawdata_L001_JFV-14_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-14_hiseq2. &   
illuminapairedend -r 170721_SND405_A_L001_JFV-15_R2.fastq 170721_SND405_A_L001_JFV-15_R1.fastq | tee rawdata_L001_JFV-15_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-15_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-16_R2.fastq 170721_SND405_A_L001_JFV-16_R1.fastq | tee rawdata_L001_JFV-16_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-16_hiseq2. &    
illuminapairedend -r 170721_SND405_A_L001_JFV-17_R2.fastq 170721_SND405_A_L001_JFV-17_R1.fastq | tee rawdata_L001_JFV-17_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-17_hiseq2. &    
illuminapairedend -r 170721_SND405_A_L001_JFV-18_R2.fastq 170721_SND405_A_L001_JFV-18_R1.fastq | tee rawdata_L001_JFV-18_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-18_hiseq2. &   
illuminapairedend -r 170721_SND405_A_L001_JFV-19_R2.fastq 170721_SND405_A_L001_JFV-19_R1.fastq | tee rawdata_L001_JFV-19_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-19_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-20_R2.fastq 170721_SND405_A_L001_JFV-20_R1.fastq | tee rawdata_L001_JFV-20_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-20_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-21_R2.fastq 170721_SND405_A_L001_JFV-21_R1.fastq | tee rawdata_L001_JFV-21_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-21_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-22_R2.fastq 170721_SND405_A_L001_JFV-22_R1.fastq | tee rawdata_L001_JFV-22_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-22_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-23_R2.fastq 170721_SND405_A_L001_JFV-23_R1.fastq | tee rawdata_L001_JFV-23_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-23_hiseq2. &
illuminapairedend -r 170721_SND405_A_L001_JFV-24_R2.fastq 170721_SND405_A_L001_JFV-24_R1.fastq | tee rawdata_L001_JFV-24_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-24_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-13_R2.fastq 170721_SND405_A_L002_JFV-13_R1.fastq | tee rawdata_L002_JFV-13_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-13_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-14_R2.fastq 170721_SND405_A_L002_JFV-14_R1.fastq | tee rawdata_L002_JFV-14_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-14_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-15_R2.fastq 170721_SND405_A_L002_JFV-15_R1.fastq | tee rawdata_L002_JFV-15_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-15_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-16_R2.fastq 170721_SND405_A_L002_JFV-16_R1.fastq | tee rawdata_L002_JFV-16_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-16_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-17_R2.fastq 170721_SND405_A_L002_JFV-17_R1.fastq | tee rawdata_L002_JFV-17_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-17_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-18_R2.fastq 170721_SND405_A_L002_JFV-18_R1.fastq | tee rawdata_L002_JFV-18_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-18_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-19_R2.fastq 170721_SND405_A_L002_JFV-19_R1.fastq | tee rawdata_L002_JFV-19_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-19_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-20_R2.fastq 170721_SND405_A_L002_JFV-20_R1.fastq | tee rawdata_L002_JFV-20_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-20_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-21_R2.fastq 170721_SND405_A_L002_JFV-21_R1.fastq | tee rawdata_L002_JFV-21_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-21_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-22_R2.fastq 170721_SND405_A_L002_JFV-22_R1.fastq | tee rawdata_L002_JFV-22_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-22_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-23_R2.fastq 170721_SND405_A_L002_JFV-23_R1.fastq | tee rawdata_L002_JFV-23_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-23_hiseq2. &
illuminapairedend -r 170721_SND405_A_L002_JFV-24_R2.fastq 170721_SND405_A_L002_JFV-24_R1.fastq | tee rawdata_L002_JFV-24_hiseq2.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-24_hiseq2. &

touch rawdata_L001_JFV-13_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-14_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-15_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-16_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-17_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-18_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-19_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-20_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-21_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-22_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-23_hiseq2.Alignement.fastq
touch rawdata_L001_JFV-24_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-13_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-14_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-15_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-16_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-17_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-18_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-19_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-20_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-21_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-22_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-23_hiseq2.Alignement.fastq
touch rawdata_L002_JFV-24_hiseq2.Alignement.fastq

ngsfilter -t DAB13.ngsfilter -u rawdata_L001_JFV-13_hiseq2.unidentified.fastq rawdata_L001_JFV-13_hiseq2.Alignement.fastq > rawdata_L001_JFV-13_hiseq2.filtered.fastq &
ngsfilter -t DAB14.ngsfilter -u rawdata_L001_JFV-14_hiseq2.unidentified.fastq rawdata_L001_JFV-14_hiseq2.Alignement.fastq > rawdata_L001_JFV-14_hiseq2.filtered.fastq &
ngsfilter -t DAB15.ngsfilter -u rawdata_L001_JFV-15_hiseq2.unidentified.fastq rawdata_L001_JFV-15_hiseq2.Alignement.fastq > rawdata_L001_JFV-15_hiseq2.filtered.fastq &
ngsfilter -t DAB16.ngsfilter -u rawdata_L001_JFV-16_hiseq2.unidentified.fastq rawdata_L001_JFV-16_hiseq2.Alignement.fastq > rawdata_L001_JFV-16_hiseq2.filtered.fastq &
ngsfilter -t DAB17.ngsfilter -u rawdata_L001_JFV-17_hiseq2.unidentified.fastq rawdata_L001_JFV-17_hiseq2.Alignement.fastq > rawdata_L001_JFV-17_hiseq2.filtered.fastq &
ngsfilter -t DAB18.ngsfilter -u rawdata_L001_JFV-18_hiseq2.unidentified.fastq rawdata_L001_JFV-18_hiseq2.Alignement.fastq > rawdata_L001_JFV-18_hiseq2.filtered.fastq &
ngsfilter -t DAB19.ngsfilter -u rawdata_L001_JFV-19_hiseq2.unidentified.fastq rawdata_L001_JFV-19_hiseq2.Alignement.fastq > rawdata_L001_JFV-19_hiseq2.filtered.fastq &
ngsfilter -t DAB20.ngsfilter -u rawdata_L001_JFV-20_hiseq2.unidentified.fastq rawdata_L001_JFV-20_hiseq2.Alignement.fastq > rawdata_L001_JFV-20_hiseq2.filtered.fastq &
ngsfilter -t DAB21.ngsfilter -u rawdata_L001_JFV-21_hiseq2.unidentified.fastq rawdata_L001_JFV-21_hiseq2.Alignement.fastq > rawdata_L001_JFV-21_hiseq2.filtered.fastq &
ngsfilter -t DAB22.ngsfilter -u rawdata_L001_JFV-22_hiseq2.unidentified.fastq rawdata_L001_JFV-22_hiseq2.Alignement.fastq > rawdata_L001_JFV-22_hiseq2.filtered.fastq &
ngsfilter -t DAB23.ngsfilter -u rawdata_L001_JFV-23_hiseq2.unidentified.fastq rawdata_L001_JFV-23_hiseq2.Alignement.fastq > rawdata_L001_JFV-23_hiseq2.filtered.fastq &
ngsfilter -t DAB24.ngsfilter -u rawdata_L001_JFV-24_hiseq2.unidentified.fastq rawdata_L001_JFV-24_hiseq2.Alignement.fastq > rawdata_L001_JFV-24_hiseq2.filtered.fastq &
ngsfilter -t DAB13.ngsfilter -u rawdata_L002_JFV-13_hiseq2.unidentified.fastq rawdata_L002_JFV-13_hiseq2.Alignement.fastq > rawdata_L002_JFV-13_hiseq2.filtered.fastq &
ngsfilter -t DAB14.ngsfilter -u rawdata_L002_JFV-14_hiseq2.unidentified.fastq rawdata_L002_JFV-14_hiseq2.Alignement.fastq > rawdata_L002_JFV-14_hiseq2.filtered.fastq &
ngsfilter -t DAB15.ngsfilter -u rawdata_L002_JFV-15_hiseq2.unidentified.fastq rawdata_L002_JFV-15_hiseq2.Alignement.fastq > rawdata_L002_JFV-15_hiseq2.filtered.fastq &
ngsfilter -t DAB16.ngsfilter -u rawdata_L002_JFV-16_hiseq2.unidentified.fastq rawdata_L002_JFV-16_hiseq2.Alignement.fastq > rawdata_L002_JFV-16_hiseq2.filtered.fastq &
ngsfilter -t DAB17.ngsfilter -u rawdata_L002_JFV-17_hiseq2.unidentified.fastq rawdata_L002_JFV-17_hiseq2.Alignement.fastq > rawdata_L002_JFV-17_hiseq2.filtered.fastq &
ngsfilter -t DAB18.ngsfilter -u rawdata_L002_JFV-18_hiseq2.unidentified.fastq rawdata_L002_JFV-18_hiseq2.Alignement.fastq > rawdata_L002_JFV-18_hiseq2.filtered.fastq &
ngsfilter -t DAB19.ngsfilter -u rawdata_L002_JFV-19_hiseq2.unidentified.fastq rawdata_L002_JFV-19_hiseq2.Alignement.fastq > rawdata_L002_JFV-19_hiseq2.filtered.fastq &
ngsfilter -t DAB20.ngsfilter -u rawdata_L002_JFV-20_hiseq2.unidentified.fastq rawdata_L002_JFV-20_hiseq2.Alignement.fastq > rawdata_L002_JFV-20_hiseq2.filtered.fastq &
ngsfilter -t DAB21.ngsfilter -u rawdata_L002_JFV-21_hiseq2.unidentified.fastq rawdata_L002_JFV-21_hiseq2.Alignement.fastq > rawdata_L002_JFV-21_hiseq2.filtered.fastq &
ngsfilter -t DAB22.ngsfilter -u rawdata_L002_JFV-22_hiseq2.unidentified.fastq rawdata_L002_JFV-22_hiseq2.Alignement.fastq > rawdata_L002_JFV-22_hiseq2.filtered.fastq &
ngsfilter -t DAB23.ngsfilter -u rawdata_L002_JFV-23_hiseq2.unidentified.fastq rawdata_L002_JFV-23_hiseq2.Alignement.fastq > rawdata_L002_JFV-23_hiseq2.filtered.fastq &
ngsfilter -t DAB24.ngsfilter -u rawdata_L002_JFV-24_hiseq2.unidentified.fastq rawdata_L002_JFV-24_hiseq2.Alignement.fastq > rawdata_L002_JFV-24_hiseq2.filtered.fastq &

cat rawdata_L001_JFV-13_hiseq2.filtered.fastq rawdata_L002_JFV-13_hiseq2.filtered.fastq > rawdata_JFV-13.filtered.fastq &
cat rawdata_L001_JFV-14_hiseq2.filtered.fastq rawdata_L002_JFV-14_hiseq2.filtered.fastq > rawdata_JFV-14.filtered.fastq &
cat rawdata_L001_JFV-15_hiseq2.filtered.fastq rawdata_L002_JFV-15_hiseq2.filtered.fastq > rawdata_JFV-15.filtered.fastq &
cat rawdata_L001_JFV-16_hiseq2.filtered.fastq rawdata_L002_JFV-16_hiseq2.filtered.fastq > rawdata_JFV-16.filtered.fastq &
cat rawdata_L001_JFV-17_hiseq2.filtered.fastq rawdata_L002_JFV-17_hiseq2.filtered.fastq > rawdata_JFV-17.filtered.fastq &
cat rawdata_L001_JFV-18_hiseq2.filtered.fastq rawdata_L002_JFV-18_hiseq2.filtered.fastq > rawdata_JFV-18.filtered.fastq &
cat rawdata_L001_JFV-19_hiseq2.filtered.fastq rawdata_L002_JFV-19_hiseq2.filtered.fastq > rawdata_JFV-19.filtered.fastq &
cat rawdata_L001_JFV-20_hiseq2.filtered.fastq rawdata_L002_JFV-20_hiseq2.filtered.fastq > rawdata_JFV-20.filtered.fastq &
cat rawdata_L001_JFV-21_hiseq2.filtered.fastq rawdata_L002_JFV-21_hiseq2.filtered.fastq > rawdata_JFV-21.filtered.fastq &
cat rawdata_L001_JFV-22_hiseq2.filtered.fastq rawdata_L002_JFV-22_hiseq2.filtered.fastq > rawdata_JFV-22.filtered.fastq &
cat rawdata_L001_JFV-23_hiseq2.filtered.fastq rawdata_L002_JFV-23_hiseq2.filtered.fastq > rawdata_JFV-23.filtered.fastq &
cat rawdata_L001_JFV-24_hiseq2.filtered.fastq rawdata_L002_JFV-24_hiseq2.filtered.fastq > rawdata_JFV-24.filtered.fastq &

obisplit -p MICROSAT.PCR_JFV-13_ -t experiment rawdata_JFV-13.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-14_ -t experiment rawdata_JFV-14.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-15_ -t experiment rawdata_JFV-15.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-16_ -t experiment rawdata_JFV-16.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-17_ -t experiment rawdata_JFV-17.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-18_ -t experiment rawdata_JFV-18.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-19_ -t experiment rawdata_JFV-19.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-20_ -t experiment rawdata_JFV-20.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-21_ -t experiment rawdata_JFV-21.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-22_ -t experiment rawdata_JFV-22.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-23_ -t experiment rawdata_JFV-23.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-24_ -t experiment rawdata_JFV-24.filtered.fastq &

obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-13_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-13_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-14_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-14_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-15_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-15_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-16_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-16_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-17_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-17_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-18_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-18_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-19_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-19_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-20_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-20_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-21_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-21_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-22_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-22_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-23_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-23_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-24_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-24_UA_MxRout1_ZF.uniq.fasta &

obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-13_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-13_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-14_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-14_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-15_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-15_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-16_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-16_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-17_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-17_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-18_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-18_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-19_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-19_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-20_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-20_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-21_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-21_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-23_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-23_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-22_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-22_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-24_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-24_UA_MxRout1_ZF.uniq.tab &





Fasteris DinAlpBear data
# download all data (5 minutes/file)
wget -r https://data.fasteris.com/private/JFV/JFV-1-12/data/ --no-parent --user=tomaz.skrbinsek@gmail.com --password=k8zS3sh1 -A fastq.gz --no-check-certificate

# move .gz files to archive in /mnt/ngs_data_storage if not already there

# This will remove all extracted .gz files, so make sure you make a copy of them before proceeding
ls *.gz | parallel gunzip

illuminapairedend -r 170714_SND405_A_L001_JFV-1_R2.fastq 170714_SND405_A_L001_JFV-1_R1.fastq | tee rawdata_L001_JFV-1_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-1_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-2_R2.fastq 170714_SND405_A_L001_JFV-2_R1.fastq | tee rawdata_L001_JFV-2_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-2_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-3_R2.fastq 170714_SND405_A_L001_JFV-3_R1.fastq | tee rawdata_L001_JFV-3_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-3_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-4_R2.fastq 170714_SND405_A_L001_JFV-4_R1.fastq | tee rawdata_L001_JFV-4_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-4_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-5_R2.fastq 170714_SND405_A_L001_JFV-5_R1.fastq | tee rawdata_L001_JFV-5_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-5_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-6_R2.fastq 170714_SND405_A_L001_JFV-6_R1.fastq | tee rawdata_L001_JFV-6_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-6_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-7_R2.fastq 170714_SND405_A_L001_JFV-7_R1.fastq | tee rawdata_L001_JFV-7_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-7_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-8_R2.fastq 170714_SND405_A_L001_JFV-8_R1.fastq | tee rawdata_L001_JFV-8_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-8_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-9_R2.fastq 170714_SND405_A_L001_JFV-9_R1.fastq | tee rawdata_L001_JFV-9_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-9_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-10_R2.fastq 170714_SND405_A_L001_JFV-10_R1.fastq | tee rawdata_L001_JFV-10_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-10_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-11_R2.fastq 170714_SND405_A_L001_JFV-11_R1.fastq | tee rawdata_L001_JFV-11_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-11_hiseq1. &
illuminapairedend -r 170714_SND405_A_L001_JFV-12_R2.fastq 170714_SND405_A_L001_JFV-12_R1.fastq | tee rawdata_L001_JFV-12_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L001_JFV-12_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-1_R2.fastq 170714_SND405_A_L002_JFV-1_R1.fastq | tee rawdata_L002_JFV-1_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-1_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-2_R2.fastq 170714_SND405_A_L002_JFV-2_R1.fastq | tee rawdata_L002_JFV-2_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-2_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-3_R2.fastq 170714_SND405_A_L002_JFV-3_R1.fastq | tee rawdata_L002_JFV-3_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-3_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-4_R2.fastq 170714_SND405_A_L002_JFV-4_R1.fastq | tee rawdata_L002_JFV-4_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-4_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-5_R2.fastq 170714_SND405_A_L002_JFV-5_R1.fastq | tee rawdata_L002_JFV-5_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-5_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-6_R2.fastq 170714_SND405_A_L002_JFV-6_R1.fastq | tee rawdata_L002_JFV-6_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-6_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-7_R2.fastq 170714_SND405_A_L002_JFV-7_R1.fastq | tee rawdata_L002_JFV-7_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-7_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-8_R2.fastq 170714_SND405_A_L002_JFV-8_R1.fastq | tee rawdata_L002_JFV-8_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-8_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-9_R2.fastq 170714_SND405_A_L002_JFV-9_R1.fastq | tee rawdata_L002_JFV-9_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-9_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-10_R2.fastq 170714_SND405_A_L002_JFV-10_R1.fastq | tee rawdata_L002_JFV-10_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-10_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-11_R2.fastq 170714_SND405_A_L002_JFV-11_R1.fastq | tee rawdata_L002_JFV-11_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-11_hiseq1. &
illuminapairedend -r 170714_SND405_A_L002_JFV-12_R2.fastq 170714_SND405_A_L002_JFV-12_R1.fastq | tee rawdata_L002_JFV-12_hiseq1.fastq | obiannotate -S goodAli:'"Alignement" if score>40.00 else "Bad"' | obisplit -t goodAli -p rawdata_L002_JFV-12_hiseq1. &

touch rawdata_L001_JFV-1_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-2_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-3_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-4_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-5_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-6_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-7_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-8_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-9_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-10_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-11_hiseq1.Alignement.fastq
touch rawdata_L001_JFV-12_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-1_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-2_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-3_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-4_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-5_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-6_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-7_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-8_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-9_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-10_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-11_hiseq1.Alignement.fastq
touch rawdata_L002_JFV-12_hiseq1.Alignement.fastq

# move or remove intermediate results (e.g. 170714_SND405_A_L002_JFV-12_R2.fastq 170714_SND405_A_L002_JFV-12_R1.fastq) if low on disk space
# feel free to run quality check on all files, a la ls *.gz | parallel -j+48 fastqc {}

ngsfilter -t DAB01.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-1.unidentified.fastq rawdata_L001_JFV-1_hiseq1.Alignement.fastq > rawdata_L001_JFV-1_hiseq1.filtered.fastq &
ngsfilter -t DAB02.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-2.unidentified.fastq rawdata_L001_JFV-2_hiseq1.Alignement.fastq > rawdata_L001_JFV-2_hiseq1.filtered.fastq &
ngsfilter -t DAB03.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-3.unidentified.fastq rawdata_L001_JFV-3_hiseq1.Alignement.fastq > rawdata_L001_JFV-3_hiseq1.filtered.fastq &
ngsfilter -t DAB04.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-4.unidentified.fastq rawdata_L001_JFV-4_hiseq1.Alignement.fastq > rawdata_L001_JFV-4_hiseq1.filtered.fastq &
ngsfilter -t DAB05.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-5.unidentified.fastq rawdata_L001_JFV-5_hiseq1.Alignement.fastq > rawdata_L001_JFV-5_hiseq1.filtered.fastq &
ngsfilter -t DAB06.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-6.unidentified.fastq rawdata_L001_JFV-6_hiseq1.Alignement.fastq > rawdata_L001_JFV-6_hiseq1.filtered.fastq &
ngsfilter -t DAB07.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-7.unidentified.fastq rawdata_L001_JFV-7_hiseq1.Alignement.fastq > rawdata_L001_JFV-7_hiseq1.filtered.fastq &
ngsfilter -t DAB08.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-8.unidentified.fastq rawdata_L001_JFV-8_hiseq1.Alignement.fastq > rawdata_L001_JFV-8_hiseq1.filtered.fastq &
ngsfilter -t DAB09.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-9.unidentified.fastq rawdata_L001_JFV-9_hiseq1.Alignement.fastq > rawdata_L001_JFV-9_hiseq1.filtered.fastq &
ngsfilter -t DAB10.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-10.unidentified.fastq rawdata_L001_JFV-10_hiseq1.Alignement.fastq > rawdata_L001_JFV-10_hiseq1.filtered.fastq &
ngsfilter -t DAB11.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-11.unidentified.fastq rawdata_L001_JFV-11_hiseq1.Alignement.fastq > rawdata_L001_JFV-11_hiseq1.filtered.fastq &
ngsfilter -t DAB12.ngsfilter -u rawdata_dab_hiseq1-L001-JFV-12.unidentified.fastq rawdata_L001_JFV-12_hiseq1.Alignement.fastq > rawdata_L001_JFV-12_hiseq1.filtered.fastq &
ngsfilter -t DAB01.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-1.unidentified.fastq rawdata_L002_JFV-1_hiseq1.Alignement.fastq > rawdata_L002_JFV-1_hiseq1.filtered.fastq &
ngsfilter -t DAB02.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-2.unidentified.fastq rawdata_L002_JFV-2_hiseq1.Alignement.fastq > rawdata_L002_JFV-2_hiseq1.filtered.fastq &
ngsfilter -t DAB03.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-3.unidentified.fastq rawdata_L002_JFV-3_hiseq1.Alignement.fastq > rawdata_L002_JFV-3_hiseq1.filtered.fastq &
ngsfilter -t DAB04.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-4.unidentified.fastq rawdata_L002_JFV-4_hiseq1.Alignement.fastq > rawdata_L002_JFV-4_hiseq1.filtered.fastq &
ngsfilter -t DAB05.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-5.unidentified.fastq rawdata_L002_JFV-5_hiseq1.Alignement.fastq > rawdata_L002_JFV-5_hiseq1.filtered.fastq &
ngsfilter -t DAB06.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-6.unidentified.fastq rawdata_L002_JFV-6_hiseq1.Alignement.fastq > rawdata_L002_JFV-6_hiseq1.filtered.fastq &
ngsfilter -t DAB07.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-7.unidentified.fastq rawdata_L002_JFV-7_hiseq1.Alignement.fastq > rawdata_L002_JFV-7_hiseq1.filtered.fastq &
ngsfilter -t DAB08.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-8.unidentified.fastq rawdata_L002_JFV-8_hiseq1.Alignement.fastq > rawdata_L002_JFV-8_hiseq1.filtered.fastq &
ngsfilter -t DAB09.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-9.unidentified.fastq rawdata_L002_JFV-9_hiseq1.Alignement.fastq > rawdata_L002_JFV-9_hiseq1.filtered.fastq &
ngsfilter -t DAB10.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-10.unidentified.fastq rawdata_L002_JFV-10_hiseq1.Alignement.fastq > rawdata_L002_JFV-10_hiseq1.filtered.fastq &
ngsfilter -t DAB11.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-11.unidentified.fastq rawdata_L002_JFV-11_hiseq1.Alignement.fastq > rawdata_L002_JFV-11_hiseq1.filtered.fastq &
ngsfilter -t DAB12.ngsfilter -u rawdata_dab_hiseq1-L002-JFV-12.unidentified.fastq rawdata_L002_JFV-12_hiseq1.Alignement.fastq > rawdata_L002_JFV-12_hiseq1.filtered.fastq &

# feel free to move/delete unidentified files
# combine runs L001 and L002
cat rawdata_L001_JFV-1_hiseq1.filtered.fastq rawdata_L002_JFV-1_hiseq1.filtered.fastq > rawdata_JFV-1.filtered.fastq &
cat rawdata_L001_JFV-2_hiseq1.filtered.fastq rawdata_L002_JFV-2_hiseq1.filtered.fastq > rawdata_JFV-2.filtered.fastq &
cat rawdata_L001_JFV-3_hiseq1.filtered.fastq rawdata_L002_JFV-3_hiseq1.filtered.fastq > rawdata_JFV-3.filtered.fastq &
cat rawdata_L001_JFV-4_hiseq1.filtered.fastq rawdata_L002_JFV-4_hiseq1.filtered.fastq > rawdata_JFV-4.filtered.fastq &
cat rawdata_L001_JFV-5_hiseq1.filtered.fastq rawdata_L002_JFV-5_hiseq1.filtered.fastq > rawdata_JFV-5.filtered.fastq &
cat rawdata_L001_JFV-6_hiseq1.filtered.fastq rawdata_L002_JFV-6_hiseq1.filtered.fastq > rawdata_JFV-6.filtered.fastq &
cat rawdata_L001_JFV-7_hiseq1.filtered.fastq rawdata_L002_JFV-7_hiseq1.filtered.fastq > rawdata_JFV-7.filtered.fastq &
cat rawdata_L001_JFV-8_hiseq1.filtered.fastq rawdata_L002_JFV-8_hiseq1.filtered.fastq > rawdata_JFV-8.filtered.fastq &
cat rawdata_L001_JFV-9_hiseq1.filtered.fastq rawdata_L002_JFV-9_hiseq1.filtered.fastq > rawdata_JFV-9.filtered.fastq &
cat rawdata_L001_JFV-10_hiseq1.filtered.fastq rawdata_L002_JFV-10_hiseq1.filtered.fastq > rawdata_JFV-10.filtered.fastq &
cat rawdata_L001_JFV-11_hiseq1.filtered.fastq rawdata_L002_JFV-11_hiseq1.filtered.fastq > rawdata_JFV-11.filtered.fastq &
cat rawdata_L001_JFV-12_hiseq1.filtered.fastq rawdata_L002_JFV-12_hiseq1.filtered.fastq > rawdata_JFV-12.filtered.fastq &

obisplit -p MICROSAT.PCR_JFV-10_ -t experiment rawdata_JFV-10.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-11_ -t experiment rawdata_JFV-11.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-12_ -t experiment rawdata_JFV-12.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-1_ -t experiment rawdata_JFV-1.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-2_ -t experiment rawdata_JFV-2.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-3_ -t experiment rawdata_JFV-3.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-4_ -t experiment rawdata_JFV-4.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-5_ -t experiment rawdata_JFV-5.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-6_ -t experiment rawdata_JFV-6.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-7_ -t experiment rawdata_JFV-7.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-8_ -t experiment rawdata_JFV-8.filtered.fastq &
obisplit -p MICROSAT.PCR_JFV-9_ -t experiment rawdata_JFV-9.filtered.fastq &

obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-10_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-10_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-11_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-11_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-12_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-12_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-1_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-1_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-2_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-2_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-3_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-3_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-4_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-4_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-5_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-5_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-6_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-6_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-7_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-7_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-8_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-8_UA_MxRout1_ZF.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_03.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_03.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_06.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_06.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_14.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_14.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_16.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_16.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_17.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_17.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_25.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_25.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_51.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_51.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_57.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_57.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_63.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_63.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_64.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_64.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_65.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_65.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_67.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_67.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_68.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_68.uniq.fasta &
obiuniq -m sample MICROSAT.PCR_JFV-9_UA_MxRout1_ZF.fastq > MICROSAT.PCR_JFV-9_UA_MxRout1_ZF.uniq.fasta &


obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-10_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-10_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-11_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-11_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-12_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-12_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-1_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-1_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-2_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-2_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-3_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-3_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-4_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-4_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-5_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-5_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-6_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-6_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-7_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-7_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-8_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-8_UA_MxRout1_ZF.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_03.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_03.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_06.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_06.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_14.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_14.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_16.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_16.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_17.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_17.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_25.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_25.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_51.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_51.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_57.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_57.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_63.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_63.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_64.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_64.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_65.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_65.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_67.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_67.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_68.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_68.uniq.tab &
obigrep -p 'count>1' MICROSAT.PCR_JFV-9_UA_MxRout1_ZF.uniq.fasta | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > MICROSAT.PCR_JFV-9_UA_MxRout1_ZF.uniq.tab &

```
