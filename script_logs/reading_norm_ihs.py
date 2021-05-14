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
data["pval_tajima"] = norm.pdf(data["standard_TajimaD"] ) #* interval

#sns.distplot(data["pval"])
sns.distplot(data["pval_tajima"])

#transform probability to z-score - both un/normalized have same z-score
data["z"] = stats.zscore(data["pval_tajima"], nan_policy="omit")
sns.distplot(data["z"])

#keep only chro, window, and z columns
data.columns
data.drop(labels=["CHROM", "BIN_START", "N_SNPS", "TajimaD", "standard_TajimaD", 
                  "pval_tajima"], axis=1, inplace = True)
data = data.set_index('window')
data.columns = ["z_D"]
data["weight_D"] = 1

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
data1["pval_fst"] = norm.pdf(data1["standard_fst"]) #* interval

#transform probability to z-score - both un/normalized have same z-score
data1["z"] = stats.zscore(data1["pval_fst"], nan_policy="omit")
sns.distplot(data1["z"])

#keep only chro, window, and z columns
data1.columns
data1.drop(labels=["CHROM", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", 
                   "MEAN_FST", "standard_fst", "pval_fst"], axis=1,
           inplace = True)
data1.columns = ["window","z_fst"]
data1 = data1.set_index('window')
data1["weighted_fst"] = 1

#merge data and data 1
result = pd.concat([data, data1], axis=1, join="outer")

###Reading iHS normalization
os.chdir(r"D:\maulana\third_project\norm_ihs\atfl")

#read file
file = "atfl_29.ihs.out.100bins.norm.10kb.windows"
data2 = pd.read_csv(file, sep="\t", header =None)
data2.describe()

#plot distribution
sns.lineplot(data2[0], data2[5])
sns.distplot(data2[5])

#standardized score
X = data2[5]
data2["standard_ihs"] = preprocessing.scale(X)

#getting probability for each point in normal distribution
data2["pval_ihs"] = norm.pdf(data2["standard_ihs"]) #* interval

#transform probability to z-score - both un/normalized have same z-score
data2["z_ihs"] = stats.zscore(data2["pval_ihs"], nan_policy="omit")
sns.distplot(data2["z_ihs"])

#keep only chro, window, and z columns
data2.columns
data2.drop(labels=[1, 2, 3, 4, 5, "standard_ihs", "pval_ihs"], axis=1,
           inplace = True)
data2.columns = ["window","z_ihs"]
data2 = data2.set_index("window")
data2["weighted_ihs"] = 1

#merge result and data 2
result = pd.concat([result, data2], axis=1, join="outer")

###Reading nSL normalization
os.chdir(r"D:\maulana\third_project\norm_nsl\atfl")

#read file
file = "atfl_29.nsl.out.100bins.norm.10kb.windows"
data3 = pd.read_csv(file, sep="\t", header =None)
data3.describe()

#plot distribution
sns.lineplot(data3[0], data3[5])
sns.distplot(data3[5])

#standardized score
X = data3[5]
data3["standard_nsl"] = preprocessing.scale(X)

#getting probability for each point in normal distribution
data3["pval_nsl"] = norm.pdf(data3["standard_nsl"]) #* interval

#transform probability to z-score - both un/normalized have same z-score
data3["z"] = stats.zscore(data3["pval_nsl"], nan_policy="omit")
sns.distplot(data3["z"])

#keep only chro, window, and z columns
data3.columns
data3.drop(labels=[1, 2, 3, 4, 5, "standard_nsl", "pval_nsl"], axis=1,
           inplace = True)
data3.columns = ["window","z_nsl"]
data3 = data3.set_index("window")
data3["weighted_nsl"] = 1

#merge result and data 3
result = pd.concat([result, data3], axis=1, join="outer")

###Reading normalization xpehh 
os.chdir(r"D:\maulana\third_project\norm_xpehh\atfl")

#read file
file = "atfl_chbi_29.xpehh.out.norm.10kb.windows"
data4 = pd.read_csv(file, sep="\t", header=None)
data4.describe()

#plot distribution
sns.lineplot(data4[0], data4[8])
sns.distplot(data4[8])

#standardized score
X = data4[8]
data4["standard_xpehh"] = preprocessing.scale(X)

#getting probability for each point in normal distribution
data4["pval_xpehh"] = norm.pdf(data4["standard_xpehh"]) #* interval

#transform probability to z-score - both un/normalized have same z-score
data4["z"] = stats.zscore(data4["pval_xpehh"], nan_policy="omit")


#keep only chro, window, and z columns
data4.columns
data4.drop(labels=[1, 2, 3, 4, 5, 6, 7, 8, "standard_xpehh", "pval_xpehh"], axis=1,
           inplace = True)
data4.columns = ["window","z_xpehh"]
data4 = data4.set_index("window")
data4["weighted_xpehh"] = 1

#merge result and data 4
result = pd.concat([result, data4], axis=1, join="outer")

##Applying meta_ss (#stop here!!)
length_col = len(result.columns)
odd = [numbers for numbers in range(length_col) if numbers % 2 == 1 ]
even = [numbers for numbers in range(length_col) if numbers % 2 == 0 ]
numerator = 0
denumerator = 0
for o, e in zip(odd, even):
    numerator = numerator + result[o] * result[e]
    denumerator = denumerator + result[e] * result[e]
    
def meta_ss (data=result):
    

def plus (x):
    if x is int:
        z = x+1
    else :
        z = print("sorry it is not formatted yet")
    return z

plus(2)
tes= [1,2,3,4,5,6,7]
plus(tes)

data = data.assign(window = lambda x: data["BIN_START"] + 1)
numbers[0] % 2 == 1:
    return [numbers[0]] + find_odds(numbers[1:])

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
