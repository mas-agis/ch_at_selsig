# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 10:59:45 2021

@author: NUWI_352019
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.getcwd()
os.chdir(r"D:\maulana\third_project/norm_ihs/atfl")
#contents of all files in the folder
dir_contents = os.listdir()
dir_contents
# initializing substring
subs = 'windows'
# using list comprehension to get string with substring 
files = [i for i in dir_contents if subs in i]
#creating empty dataframe to hold all dataset having ihs score >=4
summary = pd.DataFrame(columns=("start", "end", "ihs_score", "chr", "group"))

#reading each files and concatenate them
for file in files:
    #file = 'atfl_29.ihs.out.100bins.norm.10kb.windows'
    #slice the first four letters defining group code
    ras = file[0:4]
    #set markers for string adjacent to chr number in the file name
    char1, char2 = "_", ".ihs.out"
    #getting chromosome number by extract string char1 and char2
    chro = file[file.find(char1)+1 : file.find(char2)]
    #reading data using panda
    data = pd.read_csv(file, sep="\t", header=None, usecols=[0,1,5])
    #rename columns to match the target summary file
    data.columns = ["start", "end", "ihs_score"]
    #adding columns for chr and group
    data["chr"] = chro
    data["group"] = ras
    #append data from the current file which score >=4 to summary data
    summary = summary.append(data[data["ihs_score"]>= 4])


fig, ax = plt.subplots()
ax.plot(data["start"], data["ihs_score"], "bo" )
ax.set(xlabel='Genome', ylabel='iHS score (normalized)')
ax.grid()
plt.show()

sns.lineplot("start", "ihs_score", data=data)

data[data["ihs_score"]>= 4]
