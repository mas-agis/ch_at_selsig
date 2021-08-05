#!/bin/bash

## usage ./tabixing_grouping_vcf.sh <group code(atfl/chbt/chbi)>

group=$(echo "$1")
mapfile -t myArray < <(seq 1 29)
output="/home/masagis/data/third_project/grouped_vcf"

cd $output/$group/

for x in "${myArray[@]}"; do
	tabix -p vcf "$group"_"$x".recode.vcf.gz;
	echo "Finished chr $x at " ; date;
	done