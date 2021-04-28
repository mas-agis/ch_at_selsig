#!/bin/bash

## usage ./normalization_iHS <group code(atfl/chbt/chbi)>

group=$(echo "$1")
mapfile -t myArray < <(seq 28 29)
norm="/home/masagis/data/selscan/norm"
output="/home/masagis/data/third_project/norm_iHS"
input="/home/masagis/data/third_project/iHS"

cd $input/$group/.
files=$(ls *.out)

echo "Doing iHS normalization for $group"

for x in $files; do
	$norm --ihs --files $x --bp-win --winsize 10000 --min-snps 10 ;
	echo "Finished chr $x at " ; date; 
	done

mv *norm* $output/$group/.