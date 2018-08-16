#!/bin/bash
#Print bam file alignment statistics
#have version of samtools bc Rich is keeping 1.19 as default for Eric's pipeline for now
mkdir ../sort_flagstat
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
/usr/local/samtools-1.3/bin/samtools flagstat "$sample".bam \
> ../sort_flagstat/"$sample"_flagstat.txt \
2> ../sort_flagstat/"$sample"_flgstbam.stderr

done
