#!/bin/bash

for file in *.txt

do

sample=`echo $file | cut -f1 -d "." `

echo $sample

awk -v sample2=`basename $file| cut -f1,2,3,4,5 -d "_" ` '{print $1,$2, $4, sample2}' "$sample".pos20.cov.txt >"$sample".pos20.bed

done
