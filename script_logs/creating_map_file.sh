#!/bin/bash

## usage ./creating_map_file 

mapfile -t myArray < <(seq 29)
input="/home/masagis/data/third_project/aa_phasing"
output="/home/masagis/data/third_project/map_file"


for x in "${myArray[@]}"; do
	zcat $input/phased_"$x".gt.vcf.gz | grep -v '#' | awk '{print $1, $3, $2/1000000, $2}' > $output/map_"$x".txt
	echo "Finished chr $x at " ; date; 
	done