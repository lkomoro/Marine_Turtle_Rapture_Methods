#!/bin/bash -l

module load ngstools/20171211

bamlist=$1
out=$2
OUT_DIR=$3

FILTERS="-minQ 20 -minMapQ 20  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"
REF=/home/lk17a/scratch/GT_genome_reference_files/Cmyd.v1.1.fa
LABELS=/home/lk17a/scratch/By_sample_demultiplexed/Rapture_bam_files_NewData/CM_treelabels.short

mkdir $OUT_DIR

nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[$nInd*5/10]

#ulimit -n 1500
angsd -b $bamlist -ref $REF -out ${OUT_DIR}/$out  -GL 1 -SNP_pval 1e-12 $FILTERS -doMajorMinor 1 -doCounts 1 \
 -doMaf 2 -minMaf 0.05 -minInd $minInd -doGeno 8 -doPost 1

 NSITES=`zcat ${OUT_DIR}/$out.mafs.gz | tail -n+2 | wc -l`
 echo $NSITES

ngsDist -verbose 1 -geno ${OUT_DIR}/$out.geno.gz -probs true -n_ind $nInd -n_sites $NSITES  -o ${OUT_DIR}/$out.boot.dist -labels $LABELS -n_threads 4 -n_boot_rep 1000 -boot_block_size 500

/home/lk17a/bin/fastme/bin/fastme -D 21 -i ${OUT_DIR}/$out.boot.dist -o ${OUT_DIR}/$out.boot.tree -m b -n b

