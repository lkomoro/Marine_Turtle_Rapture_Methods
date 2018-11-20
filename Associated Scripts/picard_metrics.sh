#!/bin/bash
# percent on target reads/bases with picard test and associated metrics
#run script in directory where files are, or change path accordingly below; would need to also change cut command accordingly

module load picard/2.17.8
mkdir ../picard_metrics


for file in ./*_sortfltr.bam

do
echo $file

sample=`echo $file |cut -f1,2,3,4,5 -d "_"`

echo $sample

#input bam
BAM="$sample"_sortfltr.bam

echo $BAM

#ref genome
REF="/project/uma_lisa_komoroske/GT_genome/Cmyd.v1.1.fa"

#Bait interval list 
BLIST="/project/uma_lisa_komoroske/picard/BAITS.interval_list"

#target interval list 
TLIST="/project/uma_lisa_komoroske/picard/TARGET.interval_list"

java -jar /share/pkg/picard/2.17.8/picard.jar CollectHsMetrics I=$BAM O=../picard_metrics/"$sample"_metrics.txt R=$REF BAIT_INTERVALS=$BLIST TARGET_INTERVALS=$TLIST

done
