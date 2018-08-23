#!/bin/bash
# Convert sam file to binary format
#after aligning with bwa, moved all the .sam files to a new folder, move to this directory and run
#have version of samtools bc Rich is keeping 1.19 as default for Eric's pipeline for now
mkdir RAD_bam
for file in *.sam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample

/usr/local/samtools-1.3/bin/samtools view -Sb "$sample".sam |/usr/local/samtools-1.3/bin/samtools sort -m 16G -@6 -O bam -T temporarysort -o ./RAD_bam/"$sample"_sort.bam \
2> ./RAD_bam/"$sample"_sort.stderr

done
