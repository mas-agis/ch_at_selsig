# -*- coding: utf-8 -*-
"""
Created on Tue May 25 12:30:13 2021

@author: NUWI_352019
"""

import os
import numpy as np
import pandas as pd

#set working directory
os.getcwd()
os.chdir(r"D:\maulana\third_project\isafe\atfl")


#define home directory containing folders for isafe outputs for each breed
home_dir = r"D:\maulana\third_project\isafe"
#define target directory for binning isafe files each breed
target_dir = r"D:\maulana\third_project\norm_isafe"

class binning_isafe:
    '''Binning isafe score based on window size''' 
    #set home_dir  as base_dir
    base_dir = home_dir   
    #define class variable of chromosome length
    chr_length=(158534110, 136231102, 121005158, 120000601, 120089316, 117806340, 
            110682743, 113319770, 105454467, 103308737, 106982474, 87216183, 
            83472345, 82403003, 85007780, 81013979, 73167244, 65820629, 63449741, 
            71974595, 69862954, 60773035, 52498615, 62317253, 42350435, 51992305, 
            45612108, 45940150, 51098607) 
    #initiane instance variables of breed, tests, and chromosome       
    def __init__(self, breed, window_size):
        #intended for the name of each breed
        self.breed = str(breed)
        #intended for the size of scanning window 
        self.size = window_size
        
    def execute(self):
        for i, k in enumerate(self.chr_length):
            #define filename patterns
            filename = "final_" + str(i+1) + ".txt"
            #define the filename with its full path
            file = os.path.join(self.base_dir, self.breed, filename) 
            #read the data using pandas
            data = pd.read_csv(file, sep="\t", header= None)
            #creating bins of scanning windows
            scan_windows = np.arange(0, k, self.window_size)
            ##rename columns to match the target summary file
            data.columns = ["pos", "isafe", "daf"]
            ##(Stop in here!! pause sholat!)
            print("This is chromosome ", str(i+1), "dengan panjang ", str(k))
            
            
        #define the filename of tajima_d output using breed name and chromosome
        filename = self.breed + "_" + str(self.chr) + ".Tajima.D"


os.listdir()

#Binning outputs of iSAFE
#array of chromosome length based on ars_ucd1.2


def binning (breed="atfl", chro=29, size=10000):
    current_chr_length = chr_length[chro-1]
    scan_windows = np.arange(0, current_chr_length, size)
test = pd.read_csv("final_29.txt", sep="\t", header= None)
#rename columns to match the target summary file
    test.columns = ["pos", "isafe", "daf"]

test.groupby(pd.cut(test['pos'], bins=np.arange(0, current_chr_length, size)))['daf'].mean()
test.groupby(pd.cut(test['pos'], bins=scan_windows))['isafe', 'daf'].mean()
test
