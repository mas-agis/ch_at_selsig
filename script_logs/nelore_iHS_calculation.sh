#!/bin/bash

## usage ./iHS_calculation

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
selscan="/home/masagis/data/selscan/selscan"
map="/home/masagis/data/third_project/map_file"
output="/home/masagis/data/third_project/additional_dataset/nelore/"
input="/home/masagis/data/third_project/additional_dataset/nelore/"

echo "Calculating iHS for nelore"

cd $input

for x in "${myArray[@]}"; do
	$selscan --ihs --vcf phased_$x.gt.vcf.gz --map map_"$x".txt --out nelore_ihs-"$x" ;
	echo "Finished chr $x at " ; date; 
	done