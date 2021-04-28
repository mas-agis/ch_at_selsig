#!/bin/bash

## usage ./iHS_calculation <group code(atfl/chbt/chbi)>

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
selscan="/home/masagis/data/selscan/selscan"
output="/home/masagis/data/third_project/nSL"
input="/home/masagis/data/third_project/grouped_vcf"

echo "Calculating nSL for $group"

for x in "${myArray[@]}"; do
	$selscan --nsl --vcf $input/"$group"/"$group"_"$x".recode.vcf.gz --out $output/"$group"/"$group"_"$x" ;
	echo "Finished chr $x at " ; date; 
	done