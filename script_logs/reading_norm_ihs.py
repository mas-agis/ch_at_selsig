# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 10:59:45 2021

@author: NUWI_352019
"""
#import modules
import os
import pandas as pd
import seaborn as sns

#set working directory
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

#creating empty dataframe to hold all dataset regardless its ihs score 
all_data = pd.DataFrame(columns=("start", "end", "ihs_score", "chr", "group"))

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
    #append data from the current file to all dataset
    all_data = all_data.append(data)
    #append data from the current file which score >=4 to summary data
    summary = summary.append(data[data["ihs_score"]>= 4])

#Plotting the whole datasets 
g = sns.FacetGrid(data=all_data, col="chr", col_wrap=6, height=2)
g.map(sns.lineplot, "start", "ihs_score", alpha=0.7)


###Reading contents of Tajima's D
#import libraries
import os
import pandas as pd
import seaborn as sns

#set working directory
os.chdir(r"D:\maulana\third_project\D'_statistics\atfl")

#read Tajima's D file
file = "atfl_1.Tajima.D"
data = pd.read_csv(file, sep="\t")
data.describe()

#plot distribution
sns.lineplot("BIN_START", "TajimaD", data=data)

#confidence interval
upper_bound = data["TajimaD"].mean() + (data["TajimaD"].std()*1.5)
upper_bound
lower_bound = data["TajimaD"].mean() - (data["TajimaD"].std()*1.5)
lower_bound

#filtering data greater than upper_bound or lower than lower_bound
filtered_data = data[data["TajimaD"] > upper_bound]
filtered_data = filtered_data.append(data[data["TajimaD"] < lower_bound])

###Finding overlaps between summary(iHS) and filtered_data(Tajima's D)
result = pd.concat([df1, df4], axis=1, join="inner")

##calculating p-value for every points
import numpy as np
from scipy.stats import norm

data_start = -10
data_end = 10
data_points = 21
data = np.linspace(data_start, data_end, data_points)
data
len(data)

point_of_interest = 5
mu = np.mean(data)
sigma = np.std(data)
interval = (data_end - data_start) / (data_points - 1)

#normalized data
data_norm = (data-mu)/sigma
data_norm

#getting probability for each point in normal distribution
probability = norm.pdf(data, loc=mu, scale=sigma) #* interval
probability1 = norm.pdf(data_norm, loc=mu, scale=sigma) #* interval
probability
probability1

#assume we have have a variable with value of 12
norm.pdf(12, loc=mu, scale=sigma)
norm.pdf(22, loc=mu, scale=sigma)
help(norm)
# Specifically, norm.pdf(x, loc, scale) is identically equivalent to 
#norm.pdf(y) / scale with y = (x - loc) / scale. 

##Explanation from scipy.stats
from scipy.stats import norm
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)

mean, var, skew, kurt = norm.stats(moments='mvsk')

x = np.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)
ax.plot(x, norm.pdf(x),'r-', lw=5, alpha=0.6, label='norm pdf')

