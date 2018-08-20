#!/bin/bash -l

bamlist=$1
out=$2

nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[(1+$nInd)*5/10]

ulimit -n 1500
angsd -bam $bamlist -out $out -minQ 20 -minMapQ 10 -GL 1 -SNP_pval 1e-06 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -minInd $minInd 


