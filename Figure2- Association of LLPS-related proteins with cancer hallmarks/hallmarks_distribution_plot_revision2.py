import pandas as pd
import numpy as np
import os
import scipy.stats as stats
from scipy.stats.contingency import expected_freq
import csv
import seaborn as sns
import matplotlib.pyplot as plt

random_backgrounds_of_35_disorder_adjusted_proteins= pd.read_excel("100_35_cancer_associated_drivers_similarity_background_sets.xlsx",header=0)
# Drop the 'Set_ID' column
df_values_only = random_backgrounds_of_35_disorder_adjusted_proteins

# Flatten and remove duplicates of the uniprot IDs
random_backgrounds_of_35_disorder_adjusted_proteins_unique_ids = list(set(df_values_only.values.flatten()))
# print(len(random_backgrounds_of_35_disorder_adjusted_proteins_unique_ids))
# exit()

"""drivers"""
os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated")
cosmic = pd.read_csv("Updated_COSMIC_census_annotated.tsv", delimiter='\t', header=0)
all_cancer_genes = cosmic["UniProt_accession"].tolist()
LLPS_drivers = pd.read_excel("Final_modified_N.xlsx", sheet_name="merged(141)drivers", header=0)["uniprot"].tolist()
cancerous_drivers = list(set(LLPS_drivers) & set(all_cancer_genes))
for i, r in cosmic.iterrows():
    if str(r["UniProt_accession"]).strip() in cancerous_drivers:
        cosmic.loc[i, ["driver"]] = "driver"

os.chdir("/Users/nazanin/Desktop/LLPS_Cancer_Project_updated/hallmarks_analysis")

cosmic_determined_drivers = pd.read_excel("cosmic_determined_drivers.xlsx", sheet_name="Sheet1", header=0)
"""filtering COSMIC genes involved in the cancer hallmarks"""
cosmic_determined_drivers["COSMIC Census genes"] = 0
cosmic_determined_drivers["Cancer associated LLPS scaffolds"] = 0
cosmic_determined_drivers["disorder_adjusted_random_476_protein"] = 0

features = list(cosmic_determined_drivers.columns.values)[29:39]
# print(features)
# exit()

for index, row in cosmic_determined_drivers.iterrows():
    i = 0
    for item in features:
        if row[item] == 1 and row["UniProt_accession"].strip() in cancerous_drivers:
            i = i + 1
        cosmic_determined_drivers.loc[index, ["Cancer associated LLPS scaffolds"]] = i
for index, row in cosmic_determined_drivers.iterrows():
    j = 0
    for item in features:
        if row[item] == 1 and str(row["UniProt_accession"].strip()) not in cancerous_drivers:
            j = j + 1
        cosmic_determined_drivers.loc[index, ["COSMIC Census genes"]] = j
for index, row in cosmic_determined_drivers.iterrows():
    j = 0
    for item in features:
        if row[item] == 1 and row["UniProt_accession"].strip() in random_backgrounds_of_35_disorder_adjusted_proteins_unique_ids:
            j = j + 1
        cosmic_determined_drivers.loc[index, ["disorder_adjusted_random_476_protein"]] = j
# cosmic_determined_drivers.to_excel("data.xlsx")
# print(cosmic_determined_drivers)
# exit()

df = pd.DataFrame(columns=("Other COSMIC Census genes (684)", "Cancer-associated LLPS scaffolds (35)","Disorder matched random backgrounds (476)", 'hallmarks'))
x_labels = ["0","1", "2", "3", "4", "≥5"]
for each in x_labels:
    cancer_genes = 0
    cancer_drivers = 0
    random_BG=0
    for index, row in cosmic_determined_drivers.iterrows():
        if str(row["COSMIC Census genes"]) == each and (str(row["UniProt_accession"].strip()) not in cancerous_drivers):
            cancer_genes += 1
        if str(row["Cancer associated LLPS scaffolds"]) == each and str(row["driver"])=="driver":
            cancer_drivers += 1
        if str(row["disorder_adjusted_random_476_protein"]) == each and (str(row["UniProt_accession"].strip()) in random_backgrounds_of_35_disorder_adjusted_proteins_unique_ids):
            random_BG += 1
    if each == "≥5":
        for index, row in cosmic_determined_drivers.iterrows():
            if row["COSMIC Census genes"] == 5 or row["COSMIC Census genes"] == 6 or row["COSMIC Census genes"] == 7 or row["COSMIC Census genes"] == 8 or row[
                "COSMIC Census genes"] == 9 or row["COSMIC Census genes"] == 10:
                cancer_genes += 1
            if row["Cancer associated LLPS scaffolds"] == 5 or row["Cancer associated LLPS scaffolds"] == 6 or row["Cancer associated LLPS scaffolds"] == 7 or row["Cancer associated LLPS scaffolds"] == 8 or row[
                "Cancer associated LLPS scaffolds"] == 9 or row["Cancer associated LLPS scaffolds"] == 10:
                cancer_drivers += 1
            if (row["disorder_adjusted_random_476_protein"] == 5 or row["disorder_adjusted_random_476_protein"] == 6 or row[
                "disorder_adjusted_random_476_protein"] == 7 or row["disorder_adjusted_random_476_protein"] == 8 or row[
                "disorder_adjusted_random_476_protein"] == 9 or row["disorder_adjusted_random_476_protein"] == 10) and (str(row["UniProt_accession"].strip()) not in cancerous_drivers):
                random_BG += 1
    df.loc[each] = [cancer_genes, cancer_drivers,random_BG, each]
# print(df)
# exit()
####################################################### Normalized data
# df["COSMIC Census genes"]=(df["COSMIC Census genes"]/719)*100
# df["Cancer associated\nLLPS drivers"]=(df["Cancer associated\nLLPS drivers"]/35)*100
# df=df.set_index('hallmarks')
# for i,r in df.iterrows():
#     r["COSMIC Census genes"]=(r["COSMIC Census genes"]/53.546592)*100
#     r["Cancer associated\nLLPS drivers"]=(r["Cancer associated\nLLPS drivers"]/37.142857)*100
# print(df)
# ax = df[1:5].plot.bar(rot=0)
# ax.set_ylabel("Percentage of genes in each category")
# ax.set_xlabel("Number of assigned hallmarks")
# plt.ylim(0, 60)
# plt.savefig('figure_normalized.png')
# exit()
####################################################### raw data
df["Disorder-matched subset of COSMIC Census (476)"]=(df["Disorder matched random backgrounds (476)"]/476)*100
df["Other COSMIC Census genes (684)"]=(df["Other COSMIC Census genes (684)"]/684)*100
df["Cancer-associated LLPS scaffolds (35)"]=(df["Cancer-associated LLPS scaffolds (35)"]/35)*100

df=df.set_index('hallmarks')
# # Reorder columns: 1. COSMIC, 2. background, 3. LLPS
df = df[["Other COSMIC Census genes (684)",
         "Disorder-matched subset of COSMIC Census (476)",
         "Cancer-associated LLPS scaffolds (35)"]]

ax = df.plot.bar(rot=0,
                 color=("royalblue", "gray", "orangered"),
                 alpha=1,
                 figsize=(7, 5))

handles, labels = ax.get_legend_handles_labels()

# Desired order of labels
desired_order = [
    "Cancer-associated LLPS scaffolds (35)",
    "Other COSMIC Census genes (684)",
    "Disorder-matched subset of COSMIC Census (476)",

]

# Create a mapping from label to handle
label_to_handle = dict(zip(labels, handles))

# Reorder handles and labels
new_handles = [label_to_handle[label] for label in desired_order]
ax.legend(new_handles, desired_order, fontsize=12)
ax.set_ylabel("Percentage of genes in each category", fontsize=12)
ax.set_xlabel("Number of assigned hallmarks",fontsize=12)
ax.xaxis.set_tick_params(labelsize=12, length=5, width=2)
ax.yaxis.set_tick_params(labelsize=12, length=5, width=2)
plt.ylim(0, 100)
plt.savefig("figure_2B.png", dpi=600, bbox_inches='tight')
plt.show()
exit()
