#!/bin/bash -l

module load angsd/0.919
module load ngstools/20171211

POP_key=/home/lk17a/scratch/Rapture_scripts_and_metadata_keys/CM_POP_key
OUT_DIR="../SAF.OUT"
POP_key_list=$(cat $POP_key)
ulimit -n 1500
mkdir $OUT_DIR
FILTERS="-minQ 20 -minMapQ 20  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"
REF=/home/lk17a/scratch/GT_genome_reference_files/Cmyd.v1.1.fa #not really needed, results appear to be the same without for unfolded
ANC=/home/lk17a/scratch/GT_genome_reference_files/Cmyd.v1.1.fa #change to outgroup/ancestral state if have; (here we are just polarizing by subbing in the ref genome)

for POP in $POP_key_list
do
        echo $POP
        nInd=$(wc $POP.bamlist | awk '{print $1}')
        minInd=$[($nInd)*5/10]
        echo $minInd
        angsd -b $POP.bamlist -ref $REF -anc $ANC -out ${OUT_DIR}/${POP}.unf $FILTERS -minInd $minInd -GL 1 -doSaf 1  #&> /dev/null
        realSFS ${OUT_DIR}/${POP}.unf.saf.idx  > ${OUT_DIR}/${POP}.unf.sfs #2> /dev/null
done
