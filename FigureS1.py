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
# # Path to your txt files
# folder_path = "/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds/100_bg_samples_w_GO_sim_to_PhaSepDB_GOs"
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
#         # Make sure exactly 272 rows: truncate or pad with empty strings
#         uniprots = uniprots[:272] + [''] * (272 - len(uniprots))
#
#         columns.append(uniprots)
#
# # ✅ Create a DataFrame from the collected columns (transpose needed)
# df = pd.DataFrame(columns).transpose()
#
# # ✅ Rename columns to bg_sample_1, bg_sample_2, ...
# df.columns = [f'bg_sample_{i + 1}' for i in range(len(file_list))]
# # Save to Excel
# df.to_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds_PhaSepDB.xlsx", index_label="Row")
# exit()

####################### figure1
# First Overlap distributions based on randomized non-LLPS backgrounds
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
COSMIC=pd.read_excel("Updated_COSMIC_census_annotated.xlsx", sheet_name="Sheet1", header=0)["UniProt_accession"].dropna().tolist()
neurodegenerative_diseases=pd.read_excel("neuro_cancer.xlsx", sheet_name="Sheet1", header=0)["neuro-uniprot"].dropna().tolist()
hereditary_cancer=pd.read_excel("neuro_cancer.xlsx", sheet_name="Sheet1", header=0)["cancer-uniprot"].dropna().tolist()
other_diseases= pd.read_excel("Other diseases.xlsx", sheet_name="Sheet1", header=0)["uniprot"].dropna().tolist()


#random
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/Random_backgrounds")
df = pd.DataFrame(pd.read_excel("100_bg_samples_w_GO_sim_to_LLPSscaffolds_PhaSepDB.xlsx"))


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
## Calculating population STD
STD1 = statistics.pstdev(lst1, average1)
STD2 = statistics.pstdev(lst2, average2)
STD3 = statistics.pstdev(lst3, average3)
STD4 = statistics.pstdev(lst4, average4)
### Z-score (observed-expected)/std range
Zscore1_ = (11 - average1) / STD1
Zscore2_ = (73 - average2) / STD2
Zscore3_ = (19 - average3) / STD3
Zscore4_ = (56 - average4) / STD4
# print()
pvalue1_COSMIC = 2 * stats.norm.sf(abs(Zscore4_))
pvalue2_Neuro = 2 * stats.norm.sf(abs(Zscore3_))
pvalue3_hereditary_cancer = 2 * stats.norm.sf(abs(Zscore1_))
pvalue4_other_diseases = 2 * stats.norm.sf(abs(Zscore2_))

# Print results
print(f"Z-score hereditary_cancer based on non llps BG: {Zscore1_}, p-value: {pvalue3_hereditary_cancer}")
print(f"Z-score other_diseases: {Zscore2_}, p-value: {pvalue4_other_diseases}")
print(f"Z-score Neuro: {Zscore3_}, p-value: {pvalue2_Neuro}")
print(f"Z-score COSMIC: {Zscore4_}, p-value: {pvalue1_COSMIC}")


######################################################Second distributions
# Second Overlap distributions based on randomized disease backgrounds'
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
# fixed_uniprots = pd.read_excel("Final_modified_N.xlsx", sheet_name="Regulators (344) merged", header=0)["uniprot"].tolist()
# fixed_uniprots = pd.read_excel("Final_modified_N.xlsx", sheet_name="merged(141)drivers", header=0)["uniprot"].tolist()
fixed_uniprots = pd.read_excel("human_PhaSepDB.xlsx", sheet_name="271 unique_uniprots_human", header=0)["uniprot_entry"].tolist()
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
Zscore1_COSMIC = (56 - average1_) / STD1_
Zscore2_Neuro = (19 - average2_) / STD2_
Zscore3_hereditary_cancer = (11 - average3_) / STD3_
Zscore4_other_diseases = (73 - average4_) / STD4_
# two sided P values
pvalue1_COSMIC = 2 * stats.norm.sf(abs(Zscore1_COSMIC))
pvalue2_Neuro = 2 * stats.norm.sf(abs(Zscore2_Neuro))
pvalue3_hereditary_cancer = 2 * stats.norm.sf(abs(Zscore3_hereditary_cancer))
pvalue4_other_diseases = 2 * stats.norm.sf(abs(Zscore4_other_diseases))

# Print results
print(f"Z-score COSMIC: {Zscore1_COSMIC}, p-value: {pvalue1_COSMIC}")
print(f"Z-score Neuro: {Zscore2_Neuro}, p-value: {pvalue2_Neuro}")
print(f"Z-score Hereditary Cancer: {Zscore3_hereditary_cancer}, p-value: {pvalue3_hereditary_cancer}")
print(f"Z-score Other Diseases: {Zscore4_other_diseases}, p-value: {pvalue4_other_diseases}")



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
ax1.set_xlim(0, 80)
ax1.set_ylim(0, 300)
ax1.plot([11, 11], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [11], "Frequency": [1]})
ax1 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax1, kind='scatter', label='Overlap with true PhaSepDB scaffolds',
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
ax4.set_xlim(0,80)
ax4.set_ylim(0, 300)
ax4.plot([56, 56], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [56], "Frequency": [1]})
ax4 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax4, kind='scatter', label='Overlap with true PhaSepDB scaffolds',
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
ax2.set_xlim(0, 80)
ax2.set_ylim(0, 300)
ax2.plot([73, 73], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [73], "Frequency": [1]})
ax2 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax2, kind='scatter', label='Overlap with true PhaSepDB scaffolds',
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
ax3.set_xlim(0, 80)
ax3.set_ylim(0, 300)
ax3.plot([19, 19], [0, 50], color='red', linewidth=3)
point = pd.DataFrame({'Number of genes in the overlaps': [19], "Frequency": [1]})
ax3 = point.plot(x='Number of genes in the overlaps', y='Frequency', ax=ax3, kind='scatter', label='Overlap with true PhaSepDB scaffolds',
                 color='red', fontsize=15, marker="|", linewidth=3, s=1000)
# ax3.axvline(x = 4, color = 'r', label = 'Tested overlap',linewidth=3)
ax3.legend(loc="upper right", markerscale=0.5, scatterpoints=1, fontsize=11)
fig = plt.figure

fig.savefig("FigureS1.png")
exit()
