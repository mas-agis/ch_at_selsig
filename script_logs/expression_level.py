# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 15:27:08 2021

@author: NUWI_352019
"""

#import libraries

import os
import matplotlib.pyplot as plt
from wordcloud import WordCloud, STOPWORDS
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

#set working directory
os.chdir(r"D:\maulana\third_project\transcriptome_details_Fang_et_al")
#read FPKM - Gene Expression table
bovine_gene_atlas = pd.read_table("FPKM_Bovine_Gene_Altlas_723.txt", sep="\t")
bovine_gene_atlas.head()
bovine_gene_atlas.shape
