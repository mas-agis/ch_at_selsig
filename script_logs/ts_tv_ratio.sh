#!/bin/bash

## usage ./ts_tv_ratio 

#change directory to folder containing original un-phased vcfs
cd /home/masagis/data/third_project/snp_filter

#myArray is files contains individuals name for each breed
myArray=(re_sequences_samples_dabieshan.txt re_sequences_samples_dehong.txt re_sequences_samples_dengchuan.txt re_sequences_samples_fleckvieh.txt re_sequences_samples_fujian.txt re_sequences_samples_guanling.txt re_sequences_samples_hasake.txt re_sequences_samples_liping.txt re_sequences_samples_luxi.txt re_sequences_samples_menggu.txt re_sequences_samples_nanyang.txt re_sequences_samples_qinchuan.txt re_sequences_samples_wenling.txt re_sequences_samples_xizang.txt re_sequences_samples_yanbian.txt)

input="/home/masagis/data/third_project/snp_filter"
output="/home/masagis/data/third_project/metadata"

for x in "${myArray[@]}" 
do
	for i in {1..29}
	do 
		vcftools --gzvcf filtered_$i.vcf.gz --keep $x --TsTv-summary --out $output/$i/$x
	done
done

