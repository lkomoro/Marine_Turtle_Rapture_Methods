#!/bin/bash
# Index bam files for being able to pull out regions downstream/look at in IGV etc.
#after aligning with bwa, moved all the .sam files to a new folder, move to this directory and run
#have version of samtools bc Rich is keeping 1.19 as default for Eric's pipeline for now
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." ` \

echo $sample
/usr/local/samtools-1.3/bin/samtools index -b "$sample".bam \
> "$sample".bam.bai

done
