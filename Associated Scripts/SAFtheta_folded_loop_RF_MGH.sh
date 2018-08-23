#!/bin/bash -l

module load ngstools/20171211

POP_key=$1
OUT_FOLDER=$2
RF=$3

POP_key_list=$(cat $POP_key)
#ulimit -n 1500
mkdir $OUT_FOLDER
FILTERS="-minQ 20 -minMapQ 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"
ANC=/home/lk17a/scratch/GT_genome_reference_files/Cmyd.v1.1.fa

for POP in $POP_key_list
do
        echo $POP
        nInd=$(wc $POP.bamlist | awk '{print $1}')
        minInd=$[$nInd*5/10]
        echo $minInd
        angsd -bam $POP.bamlist -anc $ANC -out $OUT_FOLDER/${POP}.outfold $FILTERS -GL 1 -minInd $minInd -doSaf 1 -fold 1 -P 6 -rf $RF
        realSFS $OUT_FOLDER/$POP.outfold.saf.idx  > $OUT_FOLDER/$POP.outfold.sfs

  angsd -bam $POP.bamlist -out $OUT_FOLDER/$POP.outfold -doThetas 1 -doSaf 1 -pest $OUT_FOLDER/$POP.outfold.sfs -anc $ANC -GL 1 -fold 1 -rf $RF
	thetaStat do_stat $OUT_FOLDER/$POP.outfold.thetas.idx
done
