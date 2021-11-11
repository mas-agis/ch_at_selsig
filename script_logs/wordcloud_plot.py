# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 09:55:07 2021

@author: NUWI_352019
"""

#word cloud plot
import os
import matplotlib.pyplot as plt
from wordcloud import WordCloud, STOPWORDS
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

#change directory
os.getcwd()
os.chdir(r"D:\maulana\third_project")
home_dir =r"D:\maulana\third_project"

##read iHS dataset
iHS = pd.read_table("summary_unique_genes_lenovo_&_FPKM.txt", skiprows=(1),
                    usecols=(range(0,18)), na_values=("NA"))
#set column names
colname = ["GENE_FL", "FPKM_FL", "TISSUE_FL", "GENE_ME", "FPKM_ME", "TISSUE_ME",
           "GENE_YA", "FPKM_YA", "TISSUE_YA", "GENE_CHBI", "FPKM_CHBI", 
           "TISSUE_CHBI", "GENE_CHBI_M", "FPKM_CHBI_M", "TISSUE_CHBI_M", 
           "GENE_CHBI_L", "FPKM_CHBI_L", "TISSUE_CHBI_L"]
iHS.columns = colname

#word cloud plot for each FPKM
for i in range(2, 20, 3):
    comment_words = ' '
    token = [str(val) for val in iHS.iloc[:,i].dropna() ]
    tokens = [word.lower() for word in token]
    for words in tokens: 
        comment_words = comment_words + words + ' '
    wordcloud = WordCloud(min_font_size = 10, 
                          background_color ='white').generate(comment_words)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    label_title = "iHS_" + iHS.columns[i]
    plt.title(label = label_title + "(" + str(len(tokens)) + " gene(s))")
    file = os.path.join(home_dir, "fpkm", label_title)  
    plt.show()

##read iSAFE dataset
iSAFE = pd.read_table("summary_unique_genes_lenovo_&_FPKM.txt", skiprows=(1),
                    usecols=(range(18, 36)), na_values=("NA"))
#set column names
colname1 = ["GENE_FL", "FPKM_FL", "TISSUE_FL", "GENE_HA", "FPKM_HA", "TISSUE_HA",
           "GENE_ME", "FPKM_ME", "TISSUE_ME", "GENE_YA", "FPKM_YA", 
           "TISSUE_YA", "GENE_CHBI_M", "FPKM_CHBI_M", "TISSUE_CHBI_M", 
           "GENE_CHBI_L", "FPKM_CHBI_L", "TISSUE_CHBI_L"]
iSAFE.columns = colname1

#word cloud plot for each FPKM
for i in range(2, 20, 3):
    comment_words = ' '
    token = [str(val) for val in iSAFE.iloc[:,i].dropna() ]
    tokens = [word.lower() for word in token]
    for words in tokens: 
        comment_words = comment_words + words + ' '
    wordcloud = WordCloud(min_font_size = 10, 
                          background_color ='white').generate(comment_words)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    label_title = "iSAFE_" + iSAFE.columns[i]
    plt.title(label = label_title + "(" + str(len(tokens)) + " gene(s))")
    file = os.path.join(home_dir, "fpkm", label_title)  
    plt.show()
    
##read Fst dataset
Fst = pd.read_table("summary_unique_genes_lenovo_&_FPKM.txt", skiprows=(1),
                    usecols=(range(36, 48)), na_values=("NA"))
#set column names
colname2 = ["GENE_FL", "FPKM_FL", "TISSUE_FL", "GENE_HA", "FPKM_HA", "TISSUE_HA",
           "GENE_ME", "FPKM_ME", "TISSUE_ME", "GENE_YA", "FPKM_YA", 
           "TISSUE_YA"]
Fst.columns = colname2

#word cloud plot for each FPKM
for i in range(2, 13, 3):
    comment_words = ' '
    token = [str(val) for val in Fst.iloc[:,i].dropna() ]
    tokens = [word.lower() for word in token]
    for words in tokens: 
        comment_words = comment_words + words + ' '
    wordcloud = WordCloud(min_font_size = 10, 
                          background_color ='white').generate(comment_words)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    label_title = "Fst_" + Fst.columns[i]
    plt.title(label = label_title + "(" + str(len(tokens)) + " gene(s))")
    file = os.path.join(home_dir, "fpkm", label_title)  
    plt.show()
    

##plotting mean of expression
#concatenating tables of iHS, iSAFE, and Fst
result = pd.concat([iHS, iSAFE], axis=1)
result = pd.concat([result, Fst], axis=1)
#getting the mean of each FPKM column
mean_fpkm = [result.iloc[:,col].mean() for col in range(1, 49, 3)]
#label name for the FPKM column
label_mean = ["iHS_FL", "iHS_ME", "iHS_YA", "iHS_CHBI", "iHS_CHBI_M", "iHS_CHBI_L", 
              "iSAFE_FL", "iSAFE_HA", "iSAFE_ME", "iSAFE_YA", "iSAFE_CHBI_M", 
              "iSAFE_CHBI_L", "Fst_FL", "Fst_HA", "Fst_ME", "Fst_YA"]
#bar plot from the mean_fpkm 
fig, ax = plt.subplots()
ax.bar(label_mean, mean_fpkm, width=0.8, orientation = 'horizontal')
ax.set_ylabel('Mean expression')
ax.set_title('Group by breed and respective test')
ax.legend()
plt.show()
#scatter plot from the mean_fpkm 
fig, ax = plt.subplots()
ax.scatter(label_mean, mean_fpkm)
ax.set_ylabel('Mean expression')
ax.set_title('Group by breed and respective test')
ax.legend()
plt.show()

#removing outliers values 
result.iloc[75, 19] = np.nan #TFF1 FL above 10K
result.iloc[76, 19] = np.nan #TFF2 FL above 10K
result.iloc[2, 25] = np.nan #ALDOA ME above 10K
result.iloc[14, 19] = np.nan  #DES FL above 4K

#extracting values of FPKM colomns as lists
fpkm = []
for col in range(1, 49, 3):
    list1 = [ x for x in result.iloc[:,col] if str(x) != 'nan' and x < 1000]
    fpkm.append(list1)

 [x for x in countries if str(x) != 'nan']

fpkm = [result.iloc[:,col].values for col in range(1, 49, 3)]
fpkm = [(result.iloc[:,col]) for col in range(1, 49, 3)]
fpkm = pd.DataFrame(fpkm).transpose()
fpkm = float(fpkm)
#violin plot
fig, ax = plt.subplots()
ax.set_title('Default violin plot')
ax.set_ylabel('Observed FPKM')
ax.violinplot(fpkm)
ax.legend()
plt.show()

import matplotlib.pyplot as plt
import numpy as np


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')


# create test data
np.random.seed(19680801)
data = [sorted(np.random.normal(0, std, 100)) for std in range(1, 5)]

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)
data[1][10] = np.nan
ax1.set_title('Default violin plot')
ax1.set_ylabel('Observed values')
ax1.violinplot(data[1])

ax2.set_title('Customized violin plot')
parts = ax2.violinplot(
        data, showmeans=False, showmedians=False,
        showextrema=False)

for pc in parts['bodies']:
    pc.set_facecolor('#D43F3A')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

inds = np.arange(1, len(medians) + 1)
ax2.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
ax2.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
ax2.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

# set style for the axes
labels = ['A', 'B', 'C', 'D']
for ax in [ax1, ax2]:
    set_axis_style(ax, labels)

plt.subplots_adjust(bottom=0.15, wspace=0.05)
plt.show()