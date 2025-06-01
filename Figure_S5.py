import pandas as pd
import pickle
import numpy as np
import scipy.stats as stats
from scipy.stats.contingency import expected_freq
from matplotlib.ticker import FixedFormatter, FixedLocator
import csv
import statistics
import seaborn as sns
import os, glob
from matplotlib import pyplot as plt
from csv import reader
from scipy.stats import norm
import os
import pandas as pd
import glob

################################FIGURE S5
# First Overlap distributions based on randomized non-LLPS backgrounds
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
oncogenes=pd.read_excel("COSMIC_Onco_TSG.xlsx", sheet_name="Sheet1", header=0)["Onco"].dropna().tolist()
Tumor_supressor=pd.read_excel("COSMIC_Onco_TSG.xlsx", sheet_name="Sheet1", header=0)["TSG"].dropna().tolist()

#random randomized non-LLPS backgrounds
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
df = pd.DataFrame(pd.read_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds.xlsx"))

overlaps1 = []
overlaps2 = []

overlaps1 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(oncogenes) & set(randomm_89protein))
    overlaps1.append(overlap)
    i = i + 1

overlaps2 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(Tumor_supressor) & set(randomm_89protein))
    overlaps2.append(overlap)
    i = i + 1



df_ocogenes = pd.DataFrame(overlaps1, columns=["Number of genes in the overlaps"])
df_Tumor_suppressor = pd.DataFrame(overlaps2, columns=["Number of genes in the overlaps"])

#################################### Zscore calculation based on the randomized background distribution
### Retrieving columns as lists
oncogenes = df_ocogenes["Number of genes in the overlaps"].tolist()
Tumor_supressor = df_Tumor_suppressor["Number of genes in the overlaps"].tolist()

### Calculatng average
average1 = statistics.mean(oncogenes)
average2 = statistics.mean(Tumor_supressor)

### Calculating population STD
STD1 = statistics.pstdev(oncogenes, average1)
STD2 = statistics.pstdev(Tumor_supressor, average2)

### Z-score (observed-expected)/std range
Zscore1_ = (20 - average1) / STD1
Zscore2_ = (14 - average2) / STD2

Onco = 2 * stats.norm.sf(abs(Zscore1_ ))
TSG = 2 * stats.norm.sf(abs(Zscore2_))
print(Onco)
print(TSG)
print(Zscore1_)
print(Zscore2_)
# Second Overlap distributions based on randomized disease backgrounds'
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
fixed_uniprots = pd.read_excel("Final_modified_N.xlsx", sheet_name="merged(141)drivers", header=0)["uniprot"].tolist()
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
df1 = pd.DataFrame(pd.read_excel("ONCO_random_set_table.xlsx"))
df2 = pd.DataFrame(pd.read_excel("TSG_random_set_table.xlsx"))
overlaps1 = []
overlaps2 = []

j = 0
for df in [df1, df2]:
    if j == 0:
        overlaps1 = []
        i = 1
        while i < 1001:
            randomm_89protein = df1.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps1.append(overlap)
            i = i + 1
    if j == 1:
        overlaps2 = []
        i = 1
        while i < 1001:
            randomm_89protein = df2.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps2.append(overlap)
            i = i + 1
    j = j + 1

df1 = pd.DataFrame(overlaps1, columns=["Number of genes in the overlaps"])
df2 = pd.DataFrame(overlaps2, columns=["Number of genes in the overlaps"])
### Retrieving columns as lists
lst1 = df1["Number of genes in the overlaps"].tolist()
lst2 = df2["Number of genes in the overlaps"].tolist()
### Calculatng average
average1 = statistics.mean(lst1)
average2 = statistics.mean(lst2)

### Calculating population STD
STD1 = statistics.pstdev(lst1, average1)
STD2 = statistics.pstdev(lst2, average2)

### Z-score (observed-expected)/std range
Zscore1_Onco = (20 - average1) / STD1
Zscore2_TSG = (14 - average2) / STD2

pvalue1_Onco = 2 * stats.norm.sf(abs(Zscore1_Onco))
pvalue2_TSG = 2 * stats.norm.sf(abs(Zscore2_TSG))
print(pvalue1_Onco)
print(pvalue2_TSG)



################################### FIGURE
# sns.set_style("darkgrid")
# Create figure
fig = plt.figure(figsize=(13, 13))

# Create subplot axes
ax1 = fig.add_subplot(2, 2, 1)  # 1x3 grid, position 1
ax2 = fig.add_subplot(2, 2, 2)  # 1x3 grid, position 1





plt = sns.histplot(df_ocogenes["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax1)
plt = sns.histplot(df1["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax1)
ax1.lines[0].set_color('black')
ax1.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax1.set_ylabel('Frequency', fontsize=15)
ax1.set_xlim(0, 50)
ax1.set_ylim(0,350)
ax1.plot([20, 20], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [20], "Frequency": [1]})
ax1 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax1, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|", linewidth=3, s=1000)
ax1.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)


### Second plot
plt = sns.histplot(df_Tumor_suppressor["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax2)
plt = sns.histplot(df2["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax2)
ax2.lines[0].set_color('black')
ax2.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax2.set_ylabel('Frequency', fontsize=15)
ax2.set_xlim(0, 50)
ax2.set_ylim(0, 350)
ax2.plot([14, 14], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [14], "Frequency": [1]})
ax2 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax2, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|",  linewidth=3, s=1000)
# ax2.axvline(x = 11, color = 'r', label = 'Tested overlap', linewidth=3)
ax2.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)


fig = plt.figure

fig.savefig("S55.png")
exit()