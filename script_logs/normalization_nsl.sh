#!/bin/bash

## usage ./normalization_nsl <group code(atfl/chbt/chbi)>

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
norm="/home/masagis/data/selscan/norm"
output="/home/masagis/data/third_project/norm_nSL"
input="/home/masagis/data/third_project/nSL"

cd $input/$group/.
files=$(ls *.out)

echo "Doing iHS normalization for $group"

for x in $files; do
	$norm --nsl --files $x --bp-win --winsize 10000 --min-snps 10 ;
	echo "Finished chr $x at " ; date; 
	done

mv *norm* $output/$group/.