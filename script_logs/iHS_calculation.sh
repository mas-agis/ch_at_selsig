#!/bin/bash

## usage ./iHS_calculation <group code(atfl/chbt/chbi)>

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
selscan="/home/masagis/data/selscan/selscan"
map="/home/masagis/data/third_project/map_file"
output="/home/masagis/data/third_project/iHS"
input="/home/masagis/data/third_project/grouped_vcf"

echo "Calculating iHS for $group"

for x in "${myArray[@]}"; do
	$selscan --ihs --vcf $input/"$group"/"$group"_"$x".recode.vcf.gz --map $map/map_"$x".txt --out $output/"$group"/"$group"_"$x" ;
	echo "Finished chr $x at " ; date; 
	done