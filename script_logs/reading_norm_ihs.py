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

############################
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
from scipy import stats

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
probability1 = norm.pdf(data_norm, loc=0, scale=1) #* interval
probability
probability1

#transform probability to z-score - both un/normalized have same z-score
stats.zscore(probability)
stats.zscore(probability1)

###Reading contents of each statistical test output
#import libraries
import os
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import norm
from scipy import stats
from sklearn import preprocessing

###Reading contents of Tajima's D
#set working directory
os.chdir(r"D:\maulana\third_project\D'_statistics\atfl")

#read Tajima's D file
file = "atfl_29.Tajima.D"
data = pd.read_csv(file, sep="\t")
data.describe()

#create window column by adding 1 to bin-start, to match scanning window of other tests
data = data.assign(window = lambda x: data["BIN_START"] + 1)

#plot distribution
sns.lineplot("BIN_START", "TajimaD", data=data)
sns.distplot(data["TajimaD"])

#standardized score
X = data["TajimaD"]
data["standard_TajimaD"] = preprocessing.scale(X)

#getting probability for each point in normal distribution
data["pval_tajima"] = norm.pdf(data["TajimaD"]) #* interval
data["pval_std_tajima"] = norm.pdf(data["standard_TajimaD"] ) #* interval

#sns.distplot(data["pval"])
sns.distplot(data["pval_tajima"])
sns.distplot(data["pval_std_tajima"])

###Read output of Fst fst atfl_chya chr_29 
#set working directory
os.chdir(r"D:\maulana\third_project\fst\atfl")

#read file
file = "atfl_chya_29.windowed.weir.fst"
data1 = pd.read_csv(file, sep="\t")
data1.describe()

#plot distribution
sns.lineplot("BIN_START", "MEAN_FST", data=data1)
sns.distplot(data1["MEAN_FST"])

#standardized score
X = data1["MEAN_FST"]
data1["standard_fst"] = preprocessing.scale(X)

#getting probability for each point in normal distribution
data1["pval_fst"] = norm.pdf(data1["MEAN_FST"]) #* interval
data1["pval_std_fst"] = norm.pdf(data1["standard_fst"]) #* interval

###Reading iHS normalization
os.chdir(r"D:\maulana\third_project\norm_ihs\atfl")

#read file
file = "atfl_29.ihs.out.100bins.norm.10kb.windows"
data2 = pd.read_csv(file, sep="\t", header =None)
data2.describe()

#plot distribution
sns.lineplot(data2[0], data2[5])
sns.distplot(data2[5])

#getting probability for each point in normal distribution
data2["pval_ihs"] = norm.pdf(data2[5]) #* interval

###Reading nSL normalization
os.chdir(r"D:\maulana\third_project\norm_nsl\atfl")

#read file
file = "atfl_29.nsl.out.100bins.norm.10kb.windows"
data3 = pd.read_csv(file, sep="\t", header =None)
data3.describe()

#plot distribution
sns.lineplot(data3[0], data3[5])
sns.distplot(data3[5])

#getting probability for each point in normal distribution
data3["pval_nsl"] = norm.pdf(data3[5]) #* interval

###Reading normalization xpehh 
os.chdir(r"D:\maulana\third_project\norm_xpehh\atfl")

#read file
file = "atfl_chbi_29.xpehh.out.norm.10kb.windows"
data4 = pd.read_csv(file, sep="\t", header=None)
data4.describe()

#plot distribution
sns.lineplot(data4[0], data4[8])
sns.distplot(data4[8])

#getting probability for each point in normal distribution
data4["pval_xpehh"] = norm.pdf(data4[8]) #* interval

#transform probability to z-score - both un/normalized have same z-score
data["z"] = stats.zscore(data["pval"], nan_policy="omit")
data["z1"] = stats.zscore(data["pval1"], nan_policy="omit")

##making class with inputs file, which chromosomes, column of the windows/pos, 
# column of the value(default the last column of daata), 
#option whether value is p-value or not(default is yes) (5 arguments)

#data["pval"] = norm.pdf(data["TajimaD"], loc=np.mean(data["TajimaD"]), 
               #      scale=np.std(data["TajimaD"])) #* interval

#data normalization
temp = data["TajimaD"].to_numpy() 
norm = preprocessing.normalize(temp.reshape(-1, 1))

# Standardize the data attributes for the Iris dataset.
from sklearn.datasets import load_iris
from sklearn import preprocessing
# load the Iris dataset
iris = load_iris()
print(iris.data.shape)
# separate the data and target attributes
X = iris.data
y = iris.target
# standardize the data attributes
standardized_X = preprocessing.scale(X)

preprocessing.normalize(data["pval1"] )
