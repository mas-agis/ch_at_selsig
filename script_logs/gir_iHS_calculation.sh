#!/bin/bash

## usage ./iHS_calculation

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
selscan="/home/masagis/data/selscan/selscan"
map="/home/masagis/data/third_project/map_file"
output="/home/masagis/data/third_project/additional_dataset/gir/"
input="/home/masagis/data/third_project/additional_dataset/gir/"

echo "Calculating iHS for gir"

cd $input

for x in "${myArray[@]}"; do
	$selscan --ihs --vcf phased_$x.gt.vcf.gz --map map_"$x".txt --out gir_ihs-"$x" ;
	echo "Finished chr $x at " ; date; 
	done