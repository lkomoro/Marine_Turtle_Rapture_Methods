#!/bin/bash -l
module load ngstools/20171211

bamlist=$1
out=$2
OUT_DIR=$3

nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[$nInd*5/10]
mkdir $OUT_DIR

echo $nInd
echo $minInd
FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"

ulimit -n 1500
angsd -bam $bamlist -out ${OUT_DIR}/$out -minQ 20 -minMapQ 10 -GL 1 -SNP_pval 1e-06 -doMajorMinor 1 -doMaf 2 \
 -minMaf 0.05 -minInd $minInd -doGeno 32 -doPost 2 -doGlf 2 $FILTERS

gunzip ${OUT_DIR}/${out}.*.gz

count=$(sed 1d ${OUT_DIR}/${out}.mafs | wc -l | awk '{print $1}')

#ngsCovar -probfile ${OUT_DIR}/${out}.geno -outfile ${OUT_DIR}/${out}.covar -nind $nInd -nsites $count -call 0
