#!/bin/bash

## usage ./smc_chha

#change directory to folder smc
cd /home/masagis/data/third_project/smc

length=(158534110 136231102 121005158 120000601 120089316 117806340 110682743 113319770 105454467 103308737 106982474 87216183 83472345 82403003 85007780 81013979 73167244 65820629 63449741 71974595 69862954 60773035 52498615 62317253 42350435 51992305 45612108 45940150 51098607)

#Hasake 
for i in {1..29}; do docker run --rm -v $PWD:/mnt terhorst/smcpp:latest vcf2smc phased_$i.gt.vcf.gz chha/phased-$i-chha $i CHHA:Hasake_1,Hasake_2,Hasake_3,Hasake_4,Hasake_5,Hasake_6 --length ${length[$i-1]} --missing-cutoff 20000000 ; wait; docker run --rm -v $PWD:/mnt terhorst/smcpp:latest estimate 1.25e-8 --timepoints 5000 0 chha/phased-$i-chha -o chha/ ; wait; docker run --rm -v $PWD:/mnt terhorst/smcpp:latest plot chha/plot_phased-$i-chha.png chha/model.final.json -c -g 6 ;done
