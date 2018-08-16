#!/bin/bash
# run bedtools coverage for file in directory *.bam
mkdir ../bedtools_pos20_cov
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." ` \

echo $sample
bedtools coverage -hist -abam "$sample".bam -b ../RapS_SC/Rapture_sitekey_pos20NEW_2017-28-03.bed  >../bedtools_pos20_cov/"$sample".pos20.cov.txt

done
