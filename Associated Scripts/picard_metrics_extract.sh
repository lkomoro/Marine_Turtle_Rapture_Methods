#!/bin/bash
#extract the rows I need

for file in *_metrics.txt
do
#echo $file

sample=`echo $file |cut -f1,2,3,4,5 -d "_"`

#echo $sample

sed -n '8,8p;9q' "$sample"_metrics.txt  > "$sample"_metricsshort.txt 
sed -n '7,7p;7q' "$sample"_metrics.txt  > headers.txt 

done

wait

for filename in *_metricsshort.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > ./picardmetrics_combined.txt

cat headers.txt picardmetrics_combined.txt >picardmetrics_combined_wheaders.txt

awk '{print NF}' picardmetrics_combined_wheaders.txt #how many total columns
awk '{print $54 "\t" $12 "\t" $19 "\t"$20 }' picardmetrics_combined_wheaders.txt >picardmetrics_combined_wheaders_short.txt
