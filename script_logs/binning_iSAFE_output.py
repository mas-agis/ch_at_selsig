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
os.chdir(r"D:\maulana\third_project")


#define home directory containing folders for isafe outputs for each breed
home_dir = r"D:\maulana\third_project\isafe"
#define target directory for binning isafe files each breed
target_dir = r"D:\maulana\third_project\10_kb\norm_isafe"

class binning_isafe:
    '''Binning isafe score based on window size. Arguments needed is breed code
    and scanning window size ''' 
    #set home_dir  as base_dir
    base_dir = home_dir   
    #set target_dir  as out_dir
    out_dir = target_dir   
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
        self.window_size = window_size
        
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
            #rename columns to match the target summary file
            data.columns = ["pos", "isafe", "daf"]
            #binning the isafe file
            data1 = data.groupby(pd.cut(data['pos'], bins=scan_windows))['isafe', 'daf'].mean()
            #reset index and set 'pos' column as string
            data1 = data1.reset_index() 
            data1["pos"] = data1["pos"].astype(str)
            #split 'pos' column to based on characters and save it in windows dataframe
            windows = pd.DataFrame(data1.pos.str.split(r',|\(|]| ').tolist(), columns =["be","start", "sp", "end", "ed"])
            #drop unnecessary columns
            windows.drop(labels=["be", "sp", "ed", "end"], axis=1, inplace = True)
            #set windows["start"] as data1["pos"]
            data1["pos"] = windows["start"] 
            #define output filename
            filename = "bin_" + str(i+1) + "_" + str(self.window_size) 
            out_file = os.path.join(self.out_dir, self.breed, filename)
            #write csv output
            data1.to_csv(out_file, index=False)
            
#binning isafe for atfl 
atfl = binning_isafe("atfl", 10000)
atfl.execute()

#binning isafe for chbt
chbt = binning_isafe("chbt", 10000)
chbt.execute()

#binning isafe for chha
chha = binning_isafe("chha", 10000)
chha.execute()

#binning isafe for chme
chme = binning_isafe("chme", 10000)
chme.execute()

#binning isafe for chya
chya = binning_isafe("chya", 10000)
chya.execute()


##scratch
os.getcwd()
os.listdir()
os.chdir(r"D:\maulana\third_project\norm_isafe\atfl")
test = pd.read_csv("bin_29_10000", sep=",", header = "infer")            
test.isafe.describe()
test.daf.describe()
test
#create window column by adding 1 to bin-start, to match scanning window of other tests
test["pos"] = test["pos"] + 1
X = test["isafe"]
test["standard_isafe"] = preprocessing.scale(X)
#getting probability for each point in normal distribution
test["pval_isafe"] = norm.pdf(test["standard_isafe"]) 
#transform probability to z-score - both un/normalized have same z-score
test["z_isafe"] = stats.zscore(test["pval_isafe"], nan_policy="omit")
#weight of the test
test["weighted_isafe"] = 1
##ditto for daf
X = test["daf"]
test["standard_daf"] = preprocessing.scale(X)
#getting probability for each point in normal distribution
test["pval_daf"] = norm.pdf(test["standard_daf"]) 
#transform probability to z-score - both un/normalized have same z-score
test["z_daf"] = stats.zscore(test["pval_daf"], nan_policy="omit")
#weight of the test
test["weighted_daf"] = 1
#keep only chro, window, and z columns
test.columns
test = test[["pos", "z_isafe", "weighted_isafe", "z_daf", "weighted_daf"]]
len(test)
test

 #standardized score
            X = test[8]
            test["standard_xpehh"] = preprocessing.scale(X)
            #getting probability for each point in normal distribution
            test["pval_xpehh"] = norm.pdf(test["standard_xpehh"]) 
            #transform probability to z-score - both un/normalized have same z-score
            test["z"] = stats.zscore(test["pval_xpehh"], nan_policy="omit")
            #keep only chro, window, and z columns
            test.columns
            test.drop(labels=[1, 2, 3, 4, 5, 6, 7, 8, "standard_xpehh", 
                               "pval_xpehh"], axis=1, inplace = True)
            test.columns = ["window", "z_xpehh_" + str(num)]
            test = test.set_index("window")
            test["weighted_xpehh_" + str(num)] = 1/number_of_test 
            #concatenate by column data1 to the result
            data = pd.concat([data, test], axis=1, join="outer")         
            
            
#define the filename of tajima_d output using breed name and chromosome
filename = self.breed + "_" + str(self.chr) + ".Tajima.D"
#
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

test1 = test.groupby(pd.cut(test['pos'], bins=scan_windows))['isafe', 'daf'].mean()
test1 = test1.reset_index()
test1["pos"] = test1["pos"].astype(str)
pd.DataFrame(test1.pos.str.split(', ').tolist(), columns =["start", "end"])
windows = pd.DataFrame(test1.pos.str.split(r',|\(|]| ').tolist(), columns =["be","start", "sp", "end", "ed"])
windows.drop(labels=["be", "sp", "ed", "end"], axis=1, inplace = True)
test1["pos"] = windows["start"] 
