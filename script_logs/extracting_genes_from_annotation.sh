#!/bin/bash

## usage ./extracting_genes_from_annotation
## outputting unique genes from annotation files in the folder

folder="/home/masagis/data/third_project/by_snp/annotation"

cd $folder
files=$(ls)

for x in $files; do grep -oP '(?<=Gene:).*?(?=:protein_coding)' $x | sort -u > genes-$x ; done
