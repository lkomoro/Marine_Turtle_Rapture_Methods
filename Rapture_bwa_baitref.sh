#!/bin/bash
#run script in directory where files are, or change path accordingly below; would need to also change cut command accordingly

bwa index Rapture_reference.fasta

for file in ./*_RA.fastq

do
echo $file

sample=`echo $file |cut -f1,2,3,4 -d "_"`

echo $sample

bwa mem -t 10 ./Rapture_reference.fasta \
"$sample"_RA.fastq  \
>"$sample".sam 2> "$sample".stderr

done

