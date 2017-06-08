#!/bin/bash

# ====== USER PROVIDED PARAMETERS ========
N=20 # number of cores to use
ngsfilter="rawdata_scandinavia.ngsfilter"
unidentified="unidentified_samples.fastq"

R1="rawdata_scandinavia_R1.fastq"
R2="rawdata_scandinavia_R2.fastq"
# ==== END USER PROVIDED PARAMETERS =====

fastqutils split $R1 to_illumina_R1_ $N &
fastqutils split $R2 to_illumina_R2_ $N &
wait

R1=$( ls | grep -P "to_illumina_R1_" )
R2=$( ls | grep -P "to_illumina_R2_" )

# Do alignment of two reads.
for i in $( seq $N )
do
	partR1=$( ls | grep -P "^to_illumina_R1_\\.$i\\.fastq$" )
	partR2=$( ls | grep -P "^to_illumina_R2_\\.$i\\.fastq$" )
	illuminapairedend -r $partR1 $partR2 | tee to_ngsfilter_$i.fastq | obiannotate -S goodAli:'"Align" if score>40.00 else "Bad"' | obisplit -t goodAli -p to_ngsfilter_$i. &
done

# Remove intermediate aligned files.
echo $R1 | xargs rm -r
echo $R2 | xargs rm -r
ls | grep -P "\\.seq$|\\.err$" | xargs rm -r

input=$( ls | grep -P "to_ngsfilter_[0-9]+\\.Align\\.fastq" )

# For some reasons, files need to be "touched" in order to work. Perhaps
# it's closing them...
touch $unidentified

touchedbyangel=$( ls | grep -P "to_ngsfilter_"  )
for i in $touchedbyangel
do
	touch $i
done

parallel -j$N --result filtered_{} ngsfilter -t $ngsfilter -u $unidentified {} ::: $input

# Remove intermediate results.
ls | grep -P "\\.seq$|\\.err$" | xargs rm -r

ls | grep -P "filtered_.*\\.fastq$" | xargs cat > filtered_data.fastq

# remove intermediate results to conserve space
rm to_ngsfilter*

# Split reads by locus.
filtered=$( ls | grep -P "^filtered_"  )
parallel -j$N --result to_split_{} obisplit -p MS.PCR_ -t experiment {} ::: $filtered

ls | grep -P "\\.err$|\\.seq$" | xargs rm -r
ls | grep -P "^to_split_" | xargs rm -r

# Wasn't able to parallelize this, so it's just pushing tasks into the background.
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
obiuniq -m sample MS.PCR_UA_MxRout1_ZF.fastq > MS.PCR_UA_MxRout1_ZF.uniq.fasta &
wait

uniqfiles=$( ls | grep -P ".{2}+\\.uniq\\.fasta$" )

for i in $uniqfiles
do
	# this section removes the file suffix
	# https://stackoverflow.com/questions/125281/how-do-i-remove-the-file-suffix-and-path-portion-from-a-path-string-in-bash
	tab=${i%.fasta}
	tab=${tab##*/}
	obigrep -p 'count>1' $i | obiannotate -k merged_sample -k count | obiannotate --length | obisort -r -k seq_length | obitab --no-definition --output-seq > $tab.tab
done

ls | grep -P "\\.uniq\\.fasta" | xargs rm -r

echo "Done. I think."
