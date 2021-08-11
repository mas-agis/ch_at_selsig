#!/bin/bash

## usage ./smc_chme

#change directory to folder smc
cd /home/masagis/data/third_project/smc

length=(158534110 136231102 121005158 120000601 120089316 117806340 110682743 113319770 105454467 103308737 106982474 87216183 83472345 82403003 85007780 81013979 73167244 65820629 63449741 71974595 69862954 60773035 52498615 62317253 42350435 51992305 45612108 45940150 51098607)

#Menggu 
for i in {1..29}; do docker run --rm -v $PWD:/mnt terhorst/smcpp:latest vcf2smc phased_$i.gt.vcf.gz chme/phased-$i-chme $i CHME:Menggu_1,Menggu_10,Menggu_11,Menggu_12,Menggu_2,Menggu_3,Menggu_4,Menggu_5,Menggu_6,Menggu_7,Menggu_8,Menggu_9 --length ${length[$i-1]} --missing-cutoff 20000000 ; wait; docker run --rm -v $PWD:/mnt terhorst/smcpp:latest estimate 1.25e-8 --timepoints 5000 0 chme/phased-$i-chme -o chme/ ; wait; docker run --rm -v $PWD:/mnt terhorst/smcpp:latest plot chme/plot_phased-$i-chme.png chme/model.final.json -c -g 6 ;done
