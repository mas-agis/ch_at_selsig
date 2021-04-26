#!/bin/bash

## usage ./grouping_vcf.sh <filename listing individuals>

group=$(echo "$1")
mapfile -t myArray < <(seq 29)
input="/home/masagis/data/third_project/aa_phasing"
output="/home/masagis/data/third_project/grouped_vcf"

for x in "${myArray[@]}"; do
	#echo $input/phased_"$x".gt.vcf.gz >> $output/"$group"/"$group"_"$x".txt ; 
	vcftools --gzvcf $input/phased_"$x".gt.vcf.gz --keep $1 --recode --recode-INFO-all --out $output/"$group"/"$group"_"$x";
	bgzip $output/"$group"/"$group"_"$x".recode.vcf ;
	rm $output/"$group"/"$group"_"$x".recode.vcf;
	echo "Finished chr $x at " ; date;
	done