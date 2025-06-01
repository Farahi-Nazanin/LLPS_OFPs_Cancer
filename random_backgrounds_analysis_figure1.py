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

#################### Merging all the txt file backgrounds to an excel file
# Path to your txt files
# folder_path = "/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds/100_bg_samples_w_GO_sim_to_LLPSscaffolds"
#
# # ✅ Get a sorted list of all .txt files
# file_list = sorted(glob.glob(os.path.join(folder_path, "*.txt")))
#
# # ✅ Initialize list to store columns
# columns = []
#
# # ✅ Read each file
# for idx, file_path in enumerate(file_list):
#     with open(file_path, 'r') as f:
#         uniprots = [line.strip() for line in f if line.strip()]
#
#         # Make sure exactly 140 rows: truncate or pad with empty strings
#         uniprots = uniprots[:141] + [''] * (141 - len(uniprots))
#
#         columns.append(uniprots)
#
# # ✅ Create a DataFrame from the collected columns (transpose needed)
# df = pd.DataFrame(columns).transpose()
#
# # ✅ Rename columns to bg_sample_1, bg_sample_2, ...
# df.columns = [f'bg_sample_{i + 1}' for i in range(len(file_list))]
# # Save to Excel
# df.to_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds.xlsx", index_label="Row")
# exit()

####################### figure1
"""randomized overlaps---Histogram plot"""
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
COSMIC=pd.read_excel("Updated_COSMIC_census_annotated.xlsx", sheet_name="Sheet1", header=0)["UniProt_accession"].dropna().tolist()
neurodegenerative_diseases=pd.read_excel("neuro_cancer.xlsx", sheet_name="Sheet1", header=0)["neuro-uniprot"].dropna().tolist()
hereditary_cancer=pd.read_excel("neuro_cancer.xlsx", sheet_name="Sheet1", header=0)["cancer-uniprot"].dropna().tolist()
other_diseases= pd.read_excel("Other diseases.xlsx", sheet_name="Sheet1", header=0)["uniprot"].dropna().tolist()

#random
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
df = pd.DataFrame(pd.read_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds.xlsx"))


overlaps1 = []
overlaps2 = []
overlaps3 = []
overlaps4 = []


overlaps1 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(hereditary_cancer) & set(randomm_89protein))
    overlaps1.append(overlap)
    i = i + 1

overlaps2 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(other_diseases) & set(randomm_89protein))
    overlaps2.append(overlap)
    i = i + 1

overlaps3 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(neurodegenerative_diseases) & set(randomm_89protein))
    overlaps3.append(overlap)
    i = i + 1

overlaps4 = []
i = 1
while i < 101:
    randomm_89protein = df.iloc[:, i]
    overlap = len(set(COSMIC) & set(randomm_89protein))
    overlaps4.append(overlap)
    i = i + 1


df1 = pd.DataFrame(overlaps1, columns=["Number of genes in the overlaps"])
df2 = pd.DataFrame(overlaps2, columns=["Number of genes in the overlaps"])
df3 = pd.DataFrame(overlaps3, columns=["Number of genes in the overlaps"])
df4 = pd.DataFrame(overlaps4, columns=["Number of genes in the overlaps"])
#################################### Zscore calculation based on the randomized background distribution
### Retrieving columns as lists
lst1 = df1["Number of genes in the overlaps"].tolist()
lst2 = df2["Number of genes in the overlaps"].tolist()
lst3 = df3["Number of genes in the overlaps"].tolist()
lst4 = df4["Number of genes in the overlaps"].tolist()

### Calculatng average
average1 = statistics.mean(lst1)
average2 = statistics.mean(lst2)
average3 = statistics.mean(lst3)
average4 = statistics.mean(lst4)
### Calculating population STD
STD1 = statistics.pstdev(lst1, average1)
STD2 = statistics.pstdev(lst2, average2)
STD3 = statistics.pstdev(lst3, average3)
STD4 = statistics.pstdev(lst4, average4)
### Z-score (observed-expected)/std range
Zscore1_ = (4 - average1) / STD1
Zscore2_ = (35 - average2) / STD2
Zscore3_ = (11 - average3) / STD3
Zscore4_ = (35 - average4) / STD4
# print()
pvalue1_COSMIC = 2 * stats.norm.sf(abs(Zscore1_))
pvalue2_Neuro = 2 * stats.norm.sf(abs(Zscore2_))
pvalue3_hereditary_cancer = 2 * stats.norm.sf(abs(Zscore3_))
pvalue4_other_diseases = 2 * stats.norm.sf(abs(Zscore4_))


######################################################Second distributions
"""randomized overlaps---Histogram plot"""
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
# fixed_uniprots = pd.read_excel("Final_modified_N.xlsx", sheet_name="Regulators (344) merged", header=0)["uniprot"].tolist()
fixed_uniprots = pd.read_excel("Final_modified_N.xlsx", sheet_name="merged(141)drivers", header=0)["uniprot"].tolist()
# fixed_uniprots = pd.read_excel("human_PhaSepDB.xlsx", sheet_name="271 unique_uniprots_human", header=0)["uniprot_entry"].tolist()
# print(fixed_uniprots)
# exit()
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
df1_cosmic = pd.DataFrame(pd.read_excel("COSMIC_random_set_table.xlsx"))
df2_neuro= pd.DataFrame(pd.read_excel("neurodegenerative_diseases_random_set_table.xlsx"))
df3_germline= pd.DataFrame(pd.read_excel("hereditary_cancer_random_set_table.xlsx"))
df4_other= pd.DataFrame(pd.read_excel("other_diseases_random_set_table.xlsx"))


overlaps1_ = []
overlaps2_ = []
overlaps3_ = []
overlaps4_ = []

j = 0
for df in [df1_cosmic, df2_neuro, df3_germline, df4_other]:
    if j == 0:
        overlaps1_ = []
        i = 1
        while i < 1001:
            randomm_89protein = df1_cosmic.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps1_.append(overlap)
            i = i + 1
    if j == 1:
        overlaps2_ = []
        i = 1
        while i < 1001:
            randomm_89protein = df2_neuro.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps2_.append(overlap)
            i = i + 1
    if j == 2:
        overlaps3_ = []
        i = 1
        while i < 1001:
            randomm_89protein = df3_germline.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps3_.append(overlap)
            i = i + 1
    if j == 3:
        overlaps4_ = []
        i = 1
        while i < 1001:
            randomm_89protein = df4_other.iloc[:, i]
            overlap = len(set(fixed_uniprots) & set(randomm_89protein))
            overlaps4_.append(overlap)
            i = i + 1
    j = j + 1

df1_ = pd.DataFrame(overlaps1_, columns=["Number of genes in the overlaps"])
df2_ = pd.DataFrame(overlaps2_, columns=["Number of genes in the overlaps"])
df3_ = pd.DataFrame(overlaps3_, columns=["Number of genes in the overlaps"])
df4_ = pd.DataFrame(overlaps4_, columns=["Number of genes in the overlaps"])
#################################### Zscore calculation based on the randomized background distribution
### Retrieving columns as lists
lst1_ = df1_["Number of genes in the overlaps"].tolist()
lst2_ = df2_["Number of genes in the overlaps"].tolist()
lst3_ = df3_["Number of genes in the overlaps"].tolist()
lst4_ = df4_["Number of genes in the overlaps"].tolist()


### Z-score (observed-expected)/std range
Zscore1_ = (4 - average1) / STD1
Zscore2_ = (35 - average2) / STD2
Zscore3_ = (11 - average3) / STD3
Zscore4_ = (35 - average4) / STD4

### Calculatng average
average1_ = statistics.mean(lst1_)
average2_ = statistics.mean(lst2_)
average3_ = statistics.mean(lst3_)
average4_ = statistics.mean(lst4_)
### Calculating population STD
STD1_ = statistics.pstdev(lst1_, average1_)
STD2_ = statistics.pstdev(lst2_, average2_)
STD3_ = statistics.pstdev(lst3_, average3_)
STD4_ = statistics.pstdev(lst4_, average4_)
### Z-score (observed-expected)/std range
Zscore1_COSMIC = (35 - average1_) / STD1_
Zscore2_Neuro = (11 - average2_) / STD2_
Zscore3_hereditary_cancer = (4 - average3_) / STD3_
Zscore4_other_diseases = (35 - average4_) / STD4_
# Calculate p-values (two-tailed)

#################################### FIGURE
# sns.set_style("darkgrid")

# Create figure
fig = plt.figure(figsize=(13, 13))

# Create subplot axes
ax1 = fig.add_subplot(2, 2, 1)  # 1x3 grid, position 1
ax2 = fig.add_subplot(2, 2, 2)  # 1x3 grid, position 1
ax3 = fig.add_subplot(2, 2, 3)  # 1x3 grid, position 1
ax4 = fig.add_subplot(2, 2, 4)  # 1x3 grid, position 1




plt = sns.histplot(df1["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax1)
plt = sns.histplot(df3_["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax1)
ax1.lines[0].set_color('black')
ax1.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax1.set_ylabel('Frequency', fontsize=15)
ax1.set_xlim(0, 60)
ax1.set_ylim(0, 400)
ax1.plot([4, 4], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [4], "Frequency": [1]})
ax1 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax1, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|", linewidth=3, s=1000)

ax1.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)
###### other diseases section of the figure

plt = sns.histplot(df4["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax4)
plt = sns.histplot(df1_["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax4)
ax4.lines[0].set_color('black')
ax4.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax4.set_ylabel('Frequency', fontsize=15)
ax4.set_xlim(0,60)
ax4.set_ylim(0, 400)
ax4.plot([35, 35], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [35], "Frequency": [1]})
ax4 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax4, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|", linewidth=3, s=1000)
# ax4.axvline(x = 35, color = 'r', label = 'Tested overlap', linewidth=3)
ax4.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)

### Second plot
plt = sns.histplot(df2["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax2)
plt = sns.histplot(df4_["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax2)
ax2.lines[0].set_color('black')
ax2.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax2.set_ylabel('Frequency', fontsize=15)
ax2.set_xlim(0, 60)
ax2.set_ylim(0, 400)
ax2.plot([35, 35], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [35], "Frequency": [1]})
ax2 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax2, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|",  linewidth=3, s=1000)
# ax2.axvline(x = 11, color = 'r', label = 'Tested overlap', linewidth=3)
ax2.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)

### Third subplot
plt = sns.histplot(df3["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
                   label='Overlaps with randomized non-LLPS backgrounds', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax3)
plt = sns.histplot(df2_["Number of genes in the overlaps"], color='blue', alpha=0.1, kde=True,
                   label='Overlaps with randomized disease backgrounds', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
                   ax=ax3)
ax3.lines[0].set_color('black')
ax3.set_xlabel('Number of genes in the overlaps', fontsize=15)
ax3.set_ylabel('Frequency', fontsize=15)
ax3.set_xlim(0, 60)
ax3.set_ylim(0, 400)
ax3.plot([11, 11], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [11], "Frequency": [1]})
ax3 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax3, kind='scatter', label='Overlap with true LLPS scaffolds',
                 color='red', fontsize=15, marker="|", linewidth=3, s=1000)
# ax3.axvline(x = 4, color = 'r', label = 'Tested overlap',linewidth=3)
ax3.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)
fig = plt.figure

fig.savefig("Figure1.png")
exit()

# oncogenes=pd.read_excel("COSMIC_Onco_TSG.xlsx", sheet_name="Sheet1", header=0)["Onco"].dropna().tolist()
# Tumor_supressor=pd.read_excel("COSMIC_Onco_TSG.xlsx", sheet_name="Sheet1", header=0)["TSG"].dropna().tolist()
#
# #random
# os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
# df = pd.DataFrame(pd.read_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds.xlsx"))
#
#
# overlaps1 = []
# overlaps2 = []
#
#
#
# overlaps1 = []
# i = 1
# while i < 101:
#     randomm_89protein = df.iloc[:, i]
#     overlap = len(set(oncogenes) & set(randomm_89protein))
#     overlaps1.append(overlap)
#     i = i + 1
#
# overlaps2 = []
# i = 1
# while i < 101:
#     randomm_89protein = df.iloc[:, i]
#     overlap = len(set(Tumor_supressor) & set(randomm_89protein))
#     overlaps2.append(overlap)
#     i = i + 1
#
#
#
# df1 = pd.DataFrame(overlaps1, columns=["Number of genes in the overlaps"])
# df2 = pd.DataFrame(overlaps2, columns=["Number of genes in the overlaps"])
#
# #################################### Zscore calculation based on the randomized background distribution
# ### Retrieving columns as lists
# lst1 = df1["Number of genes in the overlaps"].tolist()
# lst2 = df2["Number of genes in the overlaps"].tolist()
#
#
#
# ### Calculatng average
# average1 = statistics.mean(lst1)
# average2 = statistics.mean(lst2)
#
# ### Calculating population STD
# STD1 = statistics.pstdev(lst1, average1)
# STD2 = statistics.pstdev(lst2, average2)
#
# ### Z-score (observed-expected)/std range
# Zscore1_ = (20 - average1) / STD1
# Zscore2_ = (14 - average2) / STD2
#
# # print()
# print(Zscore1_)
# print(Zscore2_)
#
#
#
# #################################### FIGURE
# # sns.set_style("darkgrid")
#
# # Create figure
# fig = plt.figure(figsize=(13, 13))
#
# # Create subplot axes
# ax1 = fig.add_subplot(2, 2, 1)  # 1x3 grid, position 1
# ax2 = fig.add_subplot(2, 2, 2)  # 1x3 grid, position 1
#
#
#
#
#
# plt = sns.histplot(df1["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
#                    label='Randomized overlaps', bins=range(0, 1000, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
#                    ax=ax1)
# ax1.lines[0].set_color('black')
# ax1.set_xlabel('Number of genes in the overlaps', fontsize=15)
# ax1.set_ylabel('Frequency', fontsize=15)
# ax1.set_xlim(0, 50)
# ax1.set_ylim(0, 60)
# # ax1.text(25, -12, f"Z = {Zscore1_:.2f}", fontsize=12, ha='center')
# point = pd.DataFrame({'Number of genes in the overlaps': [20], "Frequency": [1]})
# ax1 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax1, kind='scatter', label='Tested overlap',
#                  color='red', fontsize=15, marker="|", linewidth=3, s=1000)
# # ax1.axvline(x = 35, color = 'r', label = 'Tested overlap', linewidth=3)
# ax1.plot([20, 20], [0, 10], color='red', linewidth=3)
# ax1.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=12)
#
#
# ### Second plot
# plt = sns.histplot(df2["Number of genes in the overlaps"], color='gray', alpha=0.3, kde=True,
#                    label='Randomized overlaps', bins=range(0, 100, 1), kde_kws=dict(cut=3), edgecolor=(1, 1, 1, .4),
#                    ax=ax2)
# ax2.lines[0].set_color('black')
# ax2.set_xlabel('Number of genes in the overlaps', fontsize=15)
# ax2.set_ylabel('Frequency', fontsize=15)
# ax2.set_xlim(0, 50)
# ax2.set_ylim(0, 60)
# point = pd.DataFrame({'Number of genes in the overlaps': [14], "Frequency": [1]})
# ax2.plot([14, 14], [0, 10], color='red', linewidth=3)
# ax2 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax2, kind='scatter', label='Tested overlap',
#                  color='red', fontsize=15, marker="|",  linewidth=3, s=1000)
# # ax2.axvline(x = 11, color = 'r', label = 'Tested overlap', linewidth=3)
# ax2.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=12)
#
#
# fig = plt.figure
#
# fig.savefig("2.png")
# exit()