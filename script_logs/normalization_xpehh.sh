#!/bin/bash

## usage ./normalization_xpehh <group code(atfl/chbt/chbi)> <group code(atfl/chbt/chbi)>

group=$(echo "$1")
group1=$(echo "$2")
mapfile -t myArray < <(seq 28 29)
norm="/home/masagis/data/selscan/norm"
output="/home/masagis/data/third_project/norm_xpehh"
input="/home/masagis/data/third_project/xpehh"

cd $input/$group/.
files=$(ls *.out)

echo "Doing xpehh normalization for $group against $group1"

for x in $files; do
	$norm --xpehh --files $x --bp-win --winsize 10000 --min-snps 10 ;
	echo "Finished chr $x at " ; date; 
	done

mv *norm* $output/$group/.