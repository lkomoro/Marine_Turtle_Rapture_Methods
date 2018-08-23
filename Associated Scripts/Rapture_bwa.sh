#!/bin/bash
#run script in directory where files are, or change path accordingly below; would need to also change cut command accordingly

# bwa index Cmyd.v1.1.fa
mkdir ../samfiles
for file in ./*_RA.fastq

do
echo $file

sample=`echo $file |cut -f1,2,3,4,5 -d "_"`

echo $sample

bwa mem -t 10 /home/lkomoroske/Reference_Genomes_etc_files/GigaDB_green_turtle_reference_genome/Cmyd.v1.1.fa \
"$sample"_RA.fastq  "$sample"_RB.fastq \
>../samfiles/"$sample".sam 2> ../samfiles/"$sample".stderr

done
