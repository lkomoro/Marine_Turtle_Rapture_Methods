#!/bin/bash -l

bamlist=$1
out=$2

nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[$nInd*5/10]

echo $nInd
echo $minInd

ulimit -n 1500
angsd -bam $bamlist -out $out -minQ 20 -minMapQ 10 -GL 1 -SNP_pval 1e-06 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -minInd $minInd -doGeno 32 -doPost 2

gunzip ${out}.*.gz

#count=$(sed 1d ${out}.mafs | wc -l | awk '{print $1}')

#ngsCovar -probfile ${out}.geno -outfile ${out}.covar -nind $nInd -nsites $count -call 0
