#!/bin/bash -l

list=$1
#enter metadata file as argument 1 in command line (BestRADBarcodesCorrected_RaptureNEWSept2017_combinedmetadata.txt)
x=1
while [ $x -le 1536 ]
do

	string="sed -n ${x}p ${list}"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4, $5, $6, $7, $8}')
	set -- $var
	c1=$1 #well barcode
	c2=$2 #plate barcode
	c3=$3 #new well ID
	c4=$4 #new plate ID
	c5=$5 #species ID
	c6=$6 #turtle ID
	c7=$7 #LABID
	c8=$8 #D_ID

#run with just RA and RB files in different subfolders; run for RA, then edit accordingly and run for RB
#run with echo to check is correct, then if all looks good remove and run without
echo mv RA/SOMM190_${c2}_RA_GG${c1}TGCAGG.fastq RA_renamed_files/${c5}_${c7}_${c8}_${c4}_${c3}_RA.fastq
#SOMM190_AATCGA_RA_GGAACCGAGATGCAGG.fastq CM_DID_pH03_wA01_RA.fastq
	x=$(( $x + 1 ))

done
