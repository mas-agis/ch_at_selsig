# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:15:47 2021

@author: NUWI_352019
"""

#import libraries
import os
import pandas as pd
import re
from scipy.stats import norm
from sklearn import preprocessing
import matplotlib.pyplot as plt

#set working directory in ~/third_project/meta_ss
os.chdir(r"D:\maulana\third_project\meta_ss")

#set home dir, where test outputs are located
base_dir = r"D:\maulana\third_project"

class meta_ss:
    '''For chromosome 1-29, combining outputs of ihs, 
    nsl, xpehh, nad iSAFE tests. Outputs are z values of meta-ss with 
    corresponding chr and pos'''

    #set base_dir as home_dir  
    home_dir = base_dir 
  
    #list of default tests
    default_tests = ["ihs", "nsl", "xpehh", "isafe"]
    
    #initiane instance variables of breed, tests, and chromosome   
    def __init__(self, breed):
        self.breed = breed
        self.tests = self.default_tests      
    
    #define function for reading the ihs result
    def ihs(self): #breed = "atfl"
        ihs = pd.DataFrame({})
        for i in range(1,30):
            file_name = self.breed + "_" + str(i) + ".ihs.out"
            file = os.path.join(self.home_dir, "ihs", self.breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[1,5])
            data = data.rename(columns={1: "pos", 5: "ihs"})
            data["chr"] = i
            data["sc_ihs"] = 1
            ihs = pd.concat([ihs, data], axis=0, join="outer")
        ihs_reindex = ihs.set_index(["chr", "pos"])
        #standardized score
        X = ihs_reindex["ihs"]
        ihs_reindex["ihs"] = preprocessing.scale(X)
        return ihs_reindex.reset_index()

    #define function for reading the nsl result
    def nsl(self): #breed = "atfl"
        nsl = pd.DataFrame({})
        for i in range(1,30):
            file_name = self.breed + "_" + str(i) + ".nsl.out"
            file = os.path.join(self.home_dir, "nSL", self.breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[1,5])
            data = data.rename(columns={1: "pos", 5: "nsl"})
            data["chr"] = i
            data["sc_nsl"] = 1
            nsl = pd.concat([nsl, data], axis=0, join="outer")
        nsl_reindex = nsl.set_index(["chr", "pos"])
        X = nsl_reindex["nsl"]
        nsl_reindex["nsl"] = preprocessing.scale(X)
        return nsl_reindex.reset_index()

    #define function for reading the xpehh result 
    def xpehh(self):
        #define the path of xpehh output for current breed
        path = os.path.join(self.home_dir, "xpehh", self.breed)
        os.listdir(path)
        xpehh_comb = pd.DataFrame({})
        for i in range(1,30):
            #define searching pattern for output file of fst with same chromosomes
            pattern = r".*" + "_" + str(i) + ".xpehh.out"
            #list several fst outputs for current breed and chromosome against other breeds
            files = [os.path.join(path, f) for f in os.listdir(path) if re.match(pattern, f)]
            #keep number of fst tests for current breed, in order to weight the score
            number_of_test = len(files)
            #creating empty dataframe for concatenate all xpehh tests
            data = pd.DataFrame({})
            data.insert(loc=0, column="pos", value="int")
            data.insert(loc=1, column="chr", value="int")
            #looping over to read each comparison test results
            for num, file in enumerate(files):
                #read the file
                data4 = pd.read_csv(file, sep="\t", usecols= ["pos", "xpehh"])
                data4["chr"] = i
                data4 = data4.rename(columns={"xpehh" : "xpehh_" + str(num)})
                data = data.merge(data4, how="outer", on=["chr", "pos"] )
                data["sc_xpehh_" + str(num)] = 1/number_of_test
            data_reindex = data.set_index(["chr", "pos"])
            xpehh_comb = pd.concat([xpehh_comb, data_reindex], axis=0, join="outer")
        for tes in range(number_of_test):
            #print(tes*2)
            X = xpehh_comb.iloc[:,tes*2] 
            xpehh_comb.iloc[:,tes*2] = preprocessing.scale(X)    
        return xpehh_comb.reset_index() 
      
    #define function for reading the isafe result 
    def isafe(self): #breed = "atfl"
        isafe = pd.DataFrame({})
        for i in range(1,30):
            file_name = "final_" + str(i) + ".txt"
            file = os.path.join(self.home_dir, "isafe", self.breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[0,2])
            data = data.rename(columns={0: "pos", 2: "isafe"})
            data["chr"] = i
            data["sc_isafe"] = 1
            isafe = pd.concat([isafe, data], axis=0, join="outer")
        isafe_reindex = isafe.set_index(["chr", "pos"])
        X = isafe_reindex["isafe"]
        isafe_reindex["isafe"] = preprocessing.scale(X)
        return isafe_reindex.reset_index()
    
    #Combine results of all tests
    def combine_tests(self):
        combined = pd.DataFrame({})
        combined.insert(0, "pos", "int")      
        combined.insert(1, "chr", "int")
        combined = combined.merge(self.ihs(), how="outer", on=["chr", "pos"] )
        combined = combined.merge(self.nsl(), how="outer", on=["chr", "pos"] )
        combined = combined.merge(self.xpehh(), how="outer", on=["chr", "pos"] )
        combined = combined.merge(self.isafe(), how="outer", on=["chr", "pos"] )
        combined = combined.set_index(["chr", "pos"])
        combined.columns
        #getting the length of columns in the combined dataframe
        length_col = len(combined.columns)
        #even is the column number containing Z score for each test
        even = [numbers for numbers in range(length_col) if numbers % 2 == 0 ]
        #odd is the column number containing weight for each test
        odd = [numbers for numbers in range(length_col) if numbers % 2 == 1 ]
        even, odd
        #fill na in each weighted score test with bfill and ffill
        for i in odd:
            combined.iloc[:, i] = combined.iloc[:, i].fillna(method="ffill")
            combined.iloc[:, i] = combined.iloc[:, i].fillna(method="bfill")
        #adding columns numerator and denumerator of to combined df
        combined.insert(length_col, "numerator", 0)        
        combined.insert(length_col+1, "denumerator", 0)
        #looping over each test to update the numerator and denumerator values
        for e, o in zip(even, odd):
            print( e, o)
            #creating temp series for multiplication of Z score and its weighted test
            temp = combined.iloc[:,e].fillna(0) * combined.iloc[:,o]
            #add the new generated score to the numerator column, with filling Na as 0
            combined["numerator"] = combined["numerator"] + temp #temp.fillna(0)
            #update the denumerator by square of the max weighted value
            combined["denumerator"] = combined["denumerator"] + combined.iloc[:,o]**2
        #final meta-ss score
        combined = combined.assign(meta_ss = lambda x: combined["numerator"] / combined["denumerator"])
        #transform meta-ss to Z-score
        temp_meta = combined["meta_ss"]
        combined["meta_ss"] = preprocessing.scale(temp_meta)
        combined = combined.reset_index()
        #final score meta_ss with chr and pos
        return combined[["chr", "pos", "meta_ss"]]

#chbi
chbi = meta_ss("chbi").combine_tests()
chbi.to_csv("chbi_meta_ss.csv", sep='\t', index =False)
#atfl
atfl = meta_ss("atfl").combine_tests()
atfl.to_csv("atfl_meta_ss.csv", sep='\t', index =False)
#chbt
chbt = meta_ss("chbt").combine_tests()
chbt.to_csv("chbt_meta_ss.csv", sep='\t', index =False)
#chha
chha = meta_ss("chha").combine_tests()
chha.to_csv("chha_meta_ss.csv", sep='\t', index =False)
#chme
chme = meta_ss("chme").combine_tests()
chme.to_csv("chme_meta_ss.csv", sep='\t', index =False)
#chya
chya = meta_ss("chya").combine_tests()
chya.to_csv("chya_meta_ss.csv", sep='\t', index =False)



####SCRATCHHH
def ihs(breed): #breed = "atfl"
        ihs = pd.DataFrame({})
        for i in range(1,30):
            file_name = breed + "_" + str(i) + ".ihs.out"
            file = os.path.join(home_dir, "ihs", breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[1,5])
            data = data.rename(columns={1: "pos", 5: "ihs"})
            data["chr"] = i
            data["sc_ihs"] = 1
            ihs = pd.concat([ihs, data], axis=0, join="outer")
        ihs_reindex = ihs.set_index(["chr", "pos"])
        #standardized score
        X = ihs_reindex["ihs"]
        ihs_reindex["ihs"] = preprocessing.scale(X)
        return ihs_reindex
def nsl(breed): #breed = "atfl"
        nsl = pd.DataFrame({})
        for i in range(1,30):
            file_name = breed + "_" + str(i) + ".nsl.out"
            file = os.path.join(home_dir, "nSL", breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[1,5])
            data = data.rename(columns={1: "pos", 5: "nsl"})
            data["chr"] = i
            data["sc_nsl"] = 1
            nsl = pd.concat([nsl, data], axis=0, join="outer")
        nsl_reindex = nsl.set_index(["chr", "pos"])
        X = nsl_reindex["nsl"]
        nsl_reindex["nsl"] = preprocessing.scale(X)
        return nsl_reindex
def xpehh(breed):
        #define the path of xpehh output for current breed
        path = os.path.join(home_dir, "xpehh", breed)
        os.listdir(path)
        xpehh_comb = pd.DataFrame({})
        for i in range(1,30):
            #define searching pattern for output file of fst with same chromosomes
            pattern = r".*" + "_" + str(i) + ".xpehh.out"
            #list several fst outputs for current breed and chromosome against other breeds
            files = [os.path.join(path, f) for f in os.listdir(path) if re.match(pattern, f)]
            #keep number of fst tests for current breed, in order to weight the score
            number_of_test = len(files)
            #creating empty dataframe for concatenate all xpehh tests
            data = pd.DataFrame({})
            data.insert(loc=0, column="pos", value="int")
            data.insert(loc=1, column="chr", value="int")
            #looping over to read each comparison test results
            for num, file in enumerate(files):
                #read the file
                data4 = pd.read_csv(file, sep="\t", usecols= ["pos", "xpehh"])
                data4["chr"] = i
                data4 = data4.rename(columns={"xpehh" : "xpehh_" + str(num)})
                data = data.merge(data4, how="outer", on=["chr", "pos"] )
                data["sc_xpehh_" + str(num)] = 1/number_of_test
            data_reindex = data.set_index(["chr", "pos"])
            xpehh_comb = pd.concat([xpehh_comb, data_reindex], axis=0, join="outer")
        for tes in range(number_of_test):
            #print(tes*2)
            X = xpehh_comb.iloc[:,tes*2] 
            xpehh_comb.iloc[:,tes*2] = preprocessing.scale(X)    
        return xpehh_comb     
def isafe(breed): #breed = "atfl"
        isafe = pd.DataFrame({})
        for i in range(1,30):
            file_name = "final_" + str(i) + ".txt"
            file = os.path.join(home_dir, "isafe", breed, file_name)
            data = pd.read_csv(file, sep="\t", header = None, usecols=[0,2])
            data = data.rename(columns={0: "pos", 2: "isafe"})
            data["chr"] = i
            data["sc_isafe"] = 1
            isafe = pd.concat([isafe, data], axis=0, join="outer")
        isafe_reindex = isafe.set_index(["chr", "pos"])
        X = isafe_reindex["isafe"]
        isafe_reindex["isafe"] = preprocessing.scale(X)
        return isafe_reindex
    
atfl_ihs = ihs("atfl") 
atfl_ihs = atfl_ihs.reset_index()
atfl_nsl = nsl("atfl") 
atfl_nsl = atfl_nsl.reset_index()
atfl_xpehh = xpehh("atfl")   
atfl_xpehh = atfl_xpehh.reset_index()
atfl_isafe = isafe("atfl")
atfl_isafe = atfl_isafe.reset_index()

combined = pd.DataFrame({})
combined.insert(0, "pos", "int")      
combined.insert(1, "chr", "int")
combined = combined.merge(atfl_ihs, how="outer", on=["chr", "pos"] )
combined = combined.merge(atfl_nsl, how="outer", on=["chr", "pos"] )
combined = combined.merge(atfl_xpehh, how="outer", on=["chr", "pos"] )
combined = combined.merge(atfl_isafe, how="outer", on=["chr", "pos"] )
combined = combined.set_index(["chr", "pos"])
combined.columns
#getting the length of columns in the combined dataframe
length_col = len(combined.columns)
#even is the column number containing Z score for each test
even = [numbers for numbers in range(length_col) if numbers % 2 == 0 ]
#odd is the column number containing weight for each test
odd = [numbers for numbers in range(length_col) if numbers % 2 == 1 ]
even, odd
#fill na in each weighted score test with bfill and ffill
for i in odd:
    combined.iloc[:, i] = combined.iloc[:, i].fillna(method="ffill")
    combined.iloc[:, i] = combined.iloc[:, i].fillna(method="bfill")
 #adding columns numerator and denumerator of to combined df
combined.insert(length_col, "numerator", 0)        
combined.insert(length_col+1, "denumerator", 0)
#looping over each test to update the numerator and denumerator values
for e, o in zip(even, odd):
    print( e, o)
    #creating temp series for multiplication of Z score and its weighted test
    temp = combined.iloc[:,e].fillna(0) * combined.iloc[:,o]
    #add the new generated score to the numerator column, with filling Na as 0
    combined["numerator"] = combined["numerator"] + temp #temp.fillna(0)
    #update the denumerator by square of the max weighted value
    combined["denumerator"] = combined["denumerator"] + combined.iloc[:,o]**2
#final meta-ss score
combined = combined.assign(meta_ss = lambda x: combined["numerator"] / combined["denumerator"])
#transform meta-ss to Z-score
temp_meta = combined["meta_ss"]
combined["meta_ss"] = preprocessing.scale(temp_meta)
combined = combined.reset_index()
#write final score meta_ss
combined[["chr", "pos", "meta_ss"]]
##Seems success

combined.head().transpose()

to_fill = combined.iloc[:, 1].mode()

combined.iloc[:, 1] = combined.iloc[:, 1].fillna(method="ffill")
combined.iloc[:, 15].unique()


atfl = tests_summary("atfl")
atfl.breed
atfl.default_tests
temp = atfl.combine_tests()

chbi = tests_summary("chbi")


  
plt.show()
X = data[5]
data["standard"] = preprocessing.scale(X)
#getting probability for each point in normal distribution
data["pval"] = norm.pdf(data["standard"] ) #seems this is both tails
plt.hist(data["pval"])
plt.show()
len(data[data["pval"]<0.01])
len(data)

