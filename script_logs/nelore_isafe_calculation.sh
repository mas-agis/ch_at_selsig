#!/bin/bash

## usage ./isafe_calculation <group code(atfl/chbt/chbi)> <chromosome>

#defined basic variables and path
aa="/home/masagis/data/Cattle/ResourceBundle/outRoot/extracted_aa"
iSAFE="/home/masagis/data/iSAFE-master/src/isafe.py"
input="/home/masagis/data/third_project/additional_dataset/nelore"
output="/home/masagis/data/third_project/additional_dataset/nelore"

#since bash is 0 based, chromosome number is deducted by one
chrom=$2-1

#array of chromosome length based on ars_ucd1.2
chr_length=(158534110
136231102
121005158
120000601
120089316
117806340
110682743
113319770
105454467
103308737
106982474
87216183
83472345
82403003
85007780
81013979
73167244
65820629
63449741
71974595
69862954
60773035
52498615
62317253
42350435
51992305
45612108
45940150
51098607) 

#assign $length as pointed by $chrom
length="${chr_length[$chrom]}"

#declare an empty $start_array 
declare -a start=()

#add splitting points for each scanning window defined as 4Mb to the $start_array
for ((i=1; i<=$length; i+=4000000)); do start+=($i); done 

#add last base to the last element $start_array
start+=($length)

#getting length of $start_array minus one(since index starts at 0)
length_index=${#start[@]}-1

#do while looping (Here supposed to be the execution of iSAFE !!)
i=0
while [[ $i -lt $length_index ]]
do
  python2 $iSAFE --input $input/phased_$2.gt.vcf.gz --output $output/$i --AA $aa/aa_$2.fa --region "$2:${start[$i]}-${start[$i+1]}" --IgnoreGaps --MinRegionSize-ps 100 
  ((i++))
done

#python2.7 ~/data/iSAFE-master/src/isafe.py --input ~/data/third_project/grouped_vcf/atfl/atfl_29.recode.vcf.gz --output ~/data/third_project_simul/isafe/atfl --AA ~/data/Cattle/ResourceBundle/outRoot/aa_29.fa --region 29:1-4000000 --IgnoreGaps 

#while looping to combine all the output of test_$i
i=0
while [[ $i -lt $length_index ]]
do 
	cat $output/$i.iSAFE.out | sed '1d' >> $output/final_$2.txt
	((i++))
done

rm $output/*.iSAFE.out
