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
import matplotlib.pyplot as plt

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

##Applying meta_ss 
#getting the length of columns in the combined dataframe
length_col = len(result.columns)
#getting the length of rows(scanning windows)
length_row = len(result)
#even is the column number containing Z score for each test
even = [numbers for numbers in range(length_col) if numbers % 2 == 0 ]
#odd is the column number containing weight for each test
odd = [numbers for numbers in range(length_col) if numbers % 2 == 1 ]
#creating empty dataframe containing numerator and denumerator of meta-ss
meta_selsig = pd.DataFrame({'numerator' : []})
meta_selsig.insert(1, "denumerator", [], allow_duplicates=False)
#looping over each test to update the numerator and denumerator values
for e, o in zip(even, odd):
    #creating temp series for multiplication of Z score and its weighted test
    temp = result.iloc[:,e] * result.iloc[:,o]
    #add the new generated score to the numerator column, with filling Na as 0
    meta_selsig["numerator"] = meta_selsig.fillna(0)["numerator"] + temp.fillna(0)
    #update the denumerator by square of weighted value, the nan values in test weight is set as 1
    meta_selsig["denumerator"] = meta_selsig.fillna(0)["denumerator"] + result.fillna(1).iloc[:,o]**2
#finalizing score of meta-ss by numerator over square-root of denumerator
meta_selsig = meta_selsig.assign(score= lambda x: meta_selsig["numerator"] / np.sqrt(meta_selsig["denumerator"]))
#calculating -log(p-value) for the score for each window
meta_selsig = meta_selsig.assign(log_pval= lambda x: -1*np.log(norm.pdf(meta_selsig["score"])))
#bonferroni threshold
bonf = -1*np.log(0.05/length_row)
#plot the distribution of meta-ss score
sns.distplot(meta_selsig["score"])
#plot the distribution of -log(pval) of meta-ss score
sns.distplot(meta_selsig["log_pval"])
#lineplot of index and -log(p-value) with bonf-treshold for horizontal line
g = sns.lineplot(meta_selsig.index, meta_selsig["log_pval"], palette="pastel" )
g.axes.axhline(bonf, ls='--')
plt.show()

##making class with inputs file, which chromosomes, column of the windows/pos, 
# column of the value(default the last column of daata), 
#option whether value is p-value or not(default is yes) (5 arguments)

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

home_dir = r"D:\maulana\third_project"
chrom = 29

class tests_summary:
    '''For every chromosome, combining selsig outputs of tajima-D, fst, ihs, 
    nsl, and xpehh tests. Outputs are z values of each test and a dataframe 
    with numerator and denumerator columns for meta-ss class'''

    #set home_dir  as base_dir
    base_dir = home_dir   
  
    #list of default tests
    #default_tests = ["tajima_d"]
    default_tests = ["tajima_d", "fst", "ihs", "nsl", "xpehh"]
    
    #initiane instance variables of breed, tests, and chromosome   
    def __init__(self, breed, chro):
        self.breed = breed
        self.chr = chro
        self.tests = self.default_tests
        
    def add_tests(self, test):
        self.tests.append(test)
    
    #Reading contents of Tajima's D
    def tajima_d(self):
        #define the filename of tajima_d output using breed name and chromosome
        filename = self.breed + "_" + str(self.chr) + ".Tajima.D"
        #define the filename with its full path
        file = os.path.join(self.base_dir, "tajima_d", self.breed, filename) 
        #read the data using pandas
        data = pd.read_csv(file, sep="\t")
        #create window column by adding 1 to bin-start, to match scanning window of other tests
        data = data.assign(window = lambda x: data["BIN_START"] + 1)
        #standardized score
        X = data["TajimaD"]
        data["standard_TajimaD"] = preprocessing.scale(X)
        #getting probability for each point in normal distribution
        data["pval_tajima"] = norm.pdf(data["standard_TajimaD"] ) 
        #transform probability to z-score - both un/normalized have same z-score
        data["z"] = stats.zscore(data["pval_tajima"], nan_policy="omit")
        #keep only chro, window, and z columns
        data.drop(labels=["CHROM", "BIN_START", "N_SNPS", "TajimaD", "standard_TajimaD", 
                  "pval_tajima"], axis=1, inplace = True)
        #set window column as index
        data = data.set_index('window')
        data.columns = ["z_D"]
        data["weighted_D"] = 1
        return data
    
    #Reading contents of Fst
    def fst(self):
        #define the path of fst output for current breed
        path = os.path.join(self.base_dir, "fst", self.breed)
        #define searching pattern for output file of fst with same chromosomes
        pattern = r".*" + str(self.chr) + ".windowed.weir.fst"
        #list several fst outputs for current breed and chromosome against other breeds
        files = [os.path.join(path, f) for f in os.listdir(path) if re.match(pattern, f)]
        #keep number of fst tests for current breed, in order to weight the score
        number_of_test = len(files)
        #creating empty dataframe for concatenate all fst tests
        data = pd.DataFrame({})
        #looping over to read each comparison test results
        for num, file in enumerate(files):
            #read the file
            data1 = pd.read_csv(file, sep="\t")
            #standardized score
            X = data1["MEAN_FST"]
            data1["standard_fst"] = preprocessing.scale(X)
            #getting probability for each point in normal distribution
            data1["pval_fst"] = norm.pdf(data1["standard_fst"])
            #transform probability to z-score - both un/normalized have same z-score
            data1["z"] = stats.zscore(data1["pval_fst"], nan_policy="omit")
            #keep only chro, window, and z columns
            data1.columns
            data1.drop(labels=["CHROM", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", 
                   "MEAN_FST", "standard_fst", "pval_fst"], axis=1, inplace = True)
            data1.columns = ["window", "z_fst"+str(num)]
            data1 = data1.set_index('window')
            data1["weighted_fst_" + str(num)] = 1/number_of_test
            #concatenate by column data1 to the result
            data = pd.concat([data, data1], axis=1, join="outer")         
        return data

    #Reading contents of window-normalized ihs 
    def ihs(self):
        #define the filename of tajima_d output using breed name and chromosome
        filename = self.breed + "_" + str(self.chr) + ".ihs.out.100bins.norm.10kb.windows"
        #define the filename with its full path
        file = os.path.join(self.base_dir, "norm_ihs", self.breed, filename) 
        #read file
        data2 = pd.read_csv(file, sep="\t", header =None)
        #standardized score
        X = data2[5]
        data2["standard_ihs"] = preprocessing.scale(X)
        #getting probability for each point in normal distribution
        data2["pval_ihs"] = norm.pdf(data2["standard_ihs"])
        #transform probability to z-score - both un/normalized have same z-score
        data2["z_ihs"] = stats.zscore(data2["pval_ihs"], nan_policy="omit")
        #keep only chro, window, and z columns
        data2.columns
        data2.drop(labels=[1, 2, 3, 4, 5, "standard_ihs", "pval_ihs"], axis=1,
                   inplace = True)
        data2.columns = ["window","z_ihs"]
        data2 = data2.set_index("window")
        data2["weighted_ihs"] = 1
        return data2

    #Reading contents of window-normalized nsl 
    def nsl(self):
        #define the filename of tajima_d output using breed name and chromosome
        filename = self.breed + "_" + str(self.chr) + ".nsl.out.100bins.norm.10kb.windows"
        #define the filename with its full path
        file = os.path.join(self.base_dir, "norm_nsl", self.breed, filename) 
        #read file
        data3 = pd.read_csv(file, sep="\t", header =None)
        #standardized score
        X = data3[5]
        data3["standard_nsl"] = preprocessing.scale(X)
        #getting probability for each point in normal distribution
        data3["pval_nsl"] = norm.pdf(data3["standard_nsl"])
        #transform probability to z-score - both un/normalized have same z-score
        data3["z_nsl"] = stats.zscore(data3["pval_nsl"], nan_policy="omit")
        #keep only chro, window, and z columns
        data3.columns
        data3.drop(labels=[1, 2, 3, 4, 5, "standard_nsl", "pval_nsl"], axis=1,
                   inplace = True)
        data3.columns = ["window","z_nsl"]
        data3 = data3.set_index("window")
        data3["weighted_nsl"] = 1
        return data3
    
    #Reading contents of xpehh
    def xpehh(self):
        #define the path of xpehh output for current breed
        path = os.path.join(self.base_dir, "norm_xpehh", self.breed)
        #define searching pattern for output file of fst with same chromosomes
        pattern = r".*" + str(self.chr) + ".xpehh.out.norm.10kb.windows"
        #list several fst outputs for current breed and chromosome against other breeds
        files = [os.path.join(path, f) for f in os.listdir(path) if re.match(pattern, f)]
        #keep number of fst tests for current breed, in order to weight the score
        number_of_test = len(files)
        #creating empty dataframe for concatenate all xpehh tests
        data = pd.DataFrame({})
        #looping over to read each comparison test results
        for num, file in enumerate(files):
            #read the file
            data4 = pd.read_csv(file, sep="\t", header=None)
            #standardized score
            X = data4[8]
            data4["standard_xpehh"] = preprocessing.scale(X)
            #getting probability for each point in normal distribution
            data4["pval_xpehh"] = norm.pdf(data4["standard_xpehh"]) 
            #transform probability to z-score - both un/normalized have same z-score
            data4["z"] = stats.zscore(data4["pval_xpehh"], nan_policy="omit")
            #keep only chro, window, and z columns
            data4.columns
            data4.drop(labels=[1, 2, 3, 4, 5, 6, 7, 8, "standard_xpehh", 
                               "pval_xpehh"], axis=1, inplace = True)
            data4.columns = ["window", "z_xpehh_" + str(num)]
            data4 = data4.set_index("window")
            data4["weighted_xpehh_" + str(num)] = 1/number_of_test 
            #concatenate by column data1 to the result
            data = pd.concat([data, data4], axis=1, join="outer")         
        return data
    
    #Combine results of all tests
    def combine_tests(self):
        #column binding results of the five tests
        combined = pd.DataFrame({})
        for test in self.tests:
            if test == "tajima_d":
                combined = pd.concat([combined, self.tajima_d()], axis=1, join="outer")
            if test == "fst":
                combined = pd.concat([combined, self.fst()], axis=1, join="outer")
            if test == "ihs":
                combined = pd.concat([combined, self.ihs()], axis=1, join="outer")
            if test == "nsl":
                combined = pd.concat([combined, self.nsl()], axis=1, join="outer")    
            if test == "xpehh":
                combined = pd.concat([combined, self.xpehh()], axis=1, join="outer")
        #Applying meta_ss 
        #getting the length of columns in the combined dataframe
        length_col = len(combined.columns)
        #even is the column number containing Z score for each test
        even = [numbers for numbers in range(length_col) if numbers % 2 == 0 ]
        #odd is the column number containing weight for each test
        odd = [numbers for numbers in range(length_col) if numbers % 2 == 1 ]
        #creating empty dataframe containing numerator and denumerator of meta-ss
        summary = pd.DataFrame({'numerator' : []})
        summary.insert(1, "denumerator", [], allow_duplicates=False)
        #looping over each test to update the numerator and denumerator values
        for e, o in zip(even, odd):
            #creating temp series for multiplication of Z score and its weighted test
            temp = combined.iloc[:,e] * combined.iloc[:,o]
            #add the new generated score to the numerator column, with filling Na as 0
            summary["numerator"] = summary.fillna(0)["numerator"] + temp.fillna(0)
            #update the denumerator by square of the max weighted value
            summary["denumerator"] = summary.fillna(0)["denumerator"] + max(combined.iloc[:,o])**2
        #assign chromosome number to chr column
        summary["chr"] = self.chr
        #reset index 
        summary = summary.reset_index()
        #get list of column names
        cols = list(summary.columns)
        #re-arranging columns
        cols = cols[-1:] + cols[:-1]
        summary = summary[cols]
        return summary#, combined

class meta_ss:
    '''Binding tests_summary of all chromosomes.''' 
    chrom = list(range(1,30))
    #initiane instance variables of breed, tests, and chromosome   
    
    def __init__(self, breed):
        #intended for the name of each breed
        self.breed = str(breed)
        #intended to store concatenaded meta_score for all chromosomes
        self.scores = {}
        #intended to store bonferroni correction value
        self.bonf = []
        #intended to store scanning window size on the tests
        self.window = []
        #intended to store regions passing the threshold of bonferroni
        self.sig_regions = {}
            
    def finalize(self):
        final_df = pd.DataFrame({})
        for chr in chrom:
            temp = tests_summary(self.breed, chr).combine_tests()        
            final_df = pd.concat([final_df, temp], axis=0)
        #finalizing score of meta-ss by numerator over square-root of denumerator
        final_df = final_df.assign(score= lambda x: final_df["numerator"] / np.sqrt(final_df["denumerator"]))
        #calculating -log(p-value) for the score for each window
        final_df = final_df.assign(min_log_pval= lambda x: -1*np.log(norm.pdf(final_df["score"])))
        #store scanning window size into self
        self.window = abs(st.mode(final_df["window"]-final_df["window"].shift(-1)))
        #getting the length of rows(scanning windows) in final_df
        length_row = len(final_df)
        #set chr and window as index of final_df
        final_df = final_df.set_index(['chr', 'window'])
        #calculating bonferroni threshold and save it in instance variable
        self.bonf = -1*np.log(0.05/length_row)
        #save the final df in instance variable
        self.scores = final_df
        return final_df
    
    def filter_regions(self):
        #save regions passing the treshold
        temp = self.scores[self.scores["min_log_pval"] > self.bonf].reset_index()
        #assign end of window column
        temp = temp.assign(end= lambda x: temp["window"] + self.window)
        #re_order columns
        cols = [-1, 0, 1, 2:]
        temp = temp[cols]
        self.sig_regions = temp
        return self.sig_regions

#Running meta_ss for atfl    
atfl = meta_ss("atfl")
atfl.finalize()
atfl.scores
atfl.filter_regions()
atfl.sig_regions


#plot the distribution of meta-ss score
sns.distplot(atfl.scores["score"])
#plot the distribution of -log(pval) of meta-ss score
sns.distplot(atfl.scores["min_log_pval"])
#lineplot of index and -log(p-value) with bonf-treshold for horizontal line
temp = atfl.scores.reset_index().reset_index()
g = sns.lineplot(temp["index"], temp["min_log_pval"], hue=temp["chr"],palette="pastel")
g.axes.axhline(atfl.bonf, ls='--')
plt.show()

g = sns.scatterplot(temp["index"], temp["min_log_pval"], hue=temp["chr"],
                    legend=None ,palette="pastel")
g.axes.axhline(atfl.bonf, ls='--')
plt.show()

atfl.chr
atfl.base_dir
atfl.number
atfl.add_number(4)
atfl.add_number(5)
atfl.add_number(8)
atfl.add_number('yoho')
del atfl.number[-1]
atfl.number
atfl.calculate()
atfl.ihs()
atfl.xpehh()

chro = 30
class trial:
    number = list(range(27,chro))
    def __init__(self, name):
        self.name = name
        
    def print_out (self):
        for num in self.number:
            print ("This is file for " + self.name + " chromosome " + str(num) )
    
chbi = trial("chbi")
chbi.print_out()
    
class simple_math:
    '''Trial of making class with basic function adding one to an integer,
    float, or a list containing integer/float'''
    
    #initiane instance variables of breed, tests, and chromosome   
    def __init__(self, number):
        self.number = number
         
    def add_number(self, numb):
        self.number.append(numb)
        
    #simple function adding one as model for building class, taking int, float, and list
    def plus (self, x):
        if type(x) is int:
            return x+1
        elif type(x) is float:
            return x+1
        elif type(x) is list:
            return [i+1 for i in x] 
        else: 
            raise ValueError("Sorry your number is not defined")
    
    #execute plus function on self.number values    
    def calculate(self):
        return self.plus(self.number)

#trial of simple_math class
tes= [1,2,3,4,5,6,7]
trial = simple_math(tes)
trial.calculate()
a = 2
simple_math.calculate(2)

#simple function adding one as model for building class, taking int, float, and list
def plus (x):
    if type(x) is int:
        return x+1
    elif type(x) is float:
        return x+1
    elif type(x) is list:
        return [i+1 for i in x] 
    else: 
        raise ValueError("Sorry your number is not defined")

tes= [1,2,3,4,5,6,7]
plus(tes)
plus(23.5)
a = 2
plus(a)

data = data.assign(window = lambda x: data["BIN_START"] + 1)
numbers[0] % 2 == 1:
    return [numbers[0]] + find_odds(numbers[1:])




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
