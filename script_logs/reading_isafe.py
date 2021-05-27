# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:28:43 2021

@author: NUWI_352019
"""


#import libraries
import os
import pandas as pd
import re
import numpy as np
import seaborn as sns
from scipy.stats import norm
from scipy import stats
from sklearn import preprocessing
import matplotlib.pyplot as plt
import statistics as st

#change directory
os.chdir(r"D:\maulana\third_project\isafe\atfl")

#read each file, add chromosome column, and concatenate to empty dataframe
combined = pd.DataFrame({})
for i in range(1, 30):
    filename = "final_" + str(i) + ".txt"
    temp = pd.read_csv(filename, sep="\t", header = None, names=["pos", "isafe", "daf"])
    temp["chr"] = str(i)
    combined = pd.concat([combined, temp], axis=0, join="outer")
    print("this is turn of ", str(i))
del temp
#set chr and pos as index 
combined = combined.set_index(['chr', 'pos'])
#plot the distribution of isafe score
sns.distplot(combined["isafe"])
#plot the distribution of daf score
sns.distplot(combined["daf"])

#adding p-val of isafe
combined["pval_isafe"] = norm.pdf(combined["isafe"])
#plot the distribution of pval_isafe score
sns.distplot(combined["pval_isafe"])
combined["z_isafe"] = stats.zscore(combined["pval_isafe"], nan_policy="omit")
sns.distplot(combined["z_isafe"])
combined["z_isafe"] = stats.zscore(combined["isafe"], nan_policy="omit")
sns.distplot(combined["z_isafe"])

X = combined["isafe"]
combined["standard_isafe"] = preprocessing.scale(X)
combined["pval_isafe"] = norm.pdf(combined["standard_isafe"])
sns.distplot(combined["pval_isafe"])
combined["z_isafe"] = stats.zscore(combined["pval_isafe"], nan_policy="omit")
sns.distplot(combined["z_isafe"])
combined["z_isafe"] = stats.zscore(combined["isafe"], nan_policy="omit")
sns.distplot(combined["z_isafe"])


#read each file, add chromosome column, and concatenate to empty dataframe
combined1 = pd.DataFrame({})
for i in range(1, 30):
    filename = "final_" + str(i) + ".txt"
    temp = pd.read_csv(filename, sep="\t", header = None, names=["pos", "isafe", "daf"])
    temp["chr"] = str(i)
    #adding p-val of isafe
    X = temp["isafe"]
    temp["standard_isafe"] = preprocessing.scale(X)
    temp["pval_isafe"] = norm.pdf(temp["standard_isafe"])
    combined1 = pd.concat([combined1, temp], axis=0, join="outer")
    print("this is turn of ", str(i))
del temp
#set chr and pos as index 
combined1 = combined1.set_index(['chr', 'pos'])
stats.skew(combined["isafe"])
stats.skew(combined["daf"])
sns.distplot(combined["pval_isafe"])
