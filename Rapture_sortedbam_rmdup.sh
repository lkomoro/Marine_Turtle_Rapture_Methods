#!/bin/bash
#Remove PCR duplicates from sorted bam files
#have version of samtools bc Rich is keeping 1.19 as default for Eric's pipeline for now

mkdir ../bam_sort_fltr
for file in *_sort.bam

do

sample=`echo $file | cut -f1,2,3,4,5 -d "_" `

echo $sample

/usr/local/samtools-1.3/bin/samtools rmdup "$sample"_sort.bam ../bam_sort_fltr/"$sample"_sortfltr.bam

done
