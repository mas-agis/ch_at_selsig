#!/bin/bash

## usage ./xpehh_calculation <1st group code(atfl/chbt/chbi)> <2nd group code(atfl/chbt/chbi)>

group1=$(echo "$1")
group2=$(echo "$2")
mapfile -t myArray < <(seq 1 27)
selscan="/home/masagis/data/selscan/selscan"
map="/home/masagis/data/third_project/map_file"
output="/home/masagis/data/third_project/xpehh"
input="/home/masagis/data/third_project/grouped_vcf"

echo "Calculating xpehh for $group1 based on $group2"

for x in "${myArray[@]}"; do
	$selscan --xpehh --vcf $input/"$group1"/"$group1"_"$x".recode.vcf.gz --vcf-ref $input/"$group2"/"$group2"_"$x".recode.vcf.gz --map $map/map_"$x".txt --out $output/"$group1"/"$group1"_"$group2"_"$x" ;
	echo "Finished chr $x at " ; date; 
	done