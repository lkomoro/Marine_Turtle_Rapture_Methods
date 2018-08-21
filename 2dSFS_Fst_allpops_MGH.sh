#!/bin/bash -l
#N.B. this runs all comparisons in the matrix, which obviously is redundant for (x*y and y*x should be the same), and not useful for x*x. But can just ignore those files, or come back and add an if/else component here

module load angsd/0.919
module load ngstools/20171211

for pop1 in *.unf.saf.idx
do
  popname1=`echo $pop1 | cut -f1 -d "." `
  echo $popname1
for pop2 in *.unf.saf.idx
do
  popname2=`echo $pop2 | cut -f1 -d "." `
    echo $popname2
realSFS $pop1 $pop2 -P 5 >"$popname1"_"$popname2".ml

done
done

wait

for pop1 in *.unf.saf.idx
do
  popname1=`echo $pop1 | cut -f1 -d "." `
  echo $popname1
for pop2 in *.unf.saf.idx
do
    popname2=`echo $pop2 | cut -f1 -d "." `
      echo $popname2

realSFS fst index $pop1 $pop2 -sfs "$popname1"_"$popname2".ml -fstout "$popname1"_"$popname2"
realSFS fst stats "$popname1"_"$popname2".fst.idx >"$popname1"_"$popname2".pwFst

done
done

wait

for i in *.pwFst
do
popnames=`echo $i | cut -f1 -d "." `
   echo $popnames 
 cat "$i"
done > CM_pwFst_list
