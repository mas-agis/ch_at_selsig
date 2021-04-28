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

#reading each files and concatenate them
file = 'atfl_29.ihs.out.100bins.norm.10kb.windows'
#slice the first four letters defining group code
ras = file[0:4]
#getting chromosome number by extract string between two strings
char1, char2 = "_", ".ihs.out"
chro = file[file.find(char1)+1 : file.find(char2)]
#reading data using panda
data = pd.read_csv(file, sep="\t", header=None, usecols=[0,1,5])
data.columns = ["start", "end", "ihs_score"]
data["chr"] = chro
data["group"] = ras
data

fig, ax = plt.subplots()
ax.plot(data["start"], data["ihs_score"], "bo" )
ax.set(xlabel='Genome', ylabel='iHS score (normalized)')
ax.grid()
plt.show()

sns.lineplot("start", "ihs_score", data=data)

data[data["ihs_score"]>= 4]
