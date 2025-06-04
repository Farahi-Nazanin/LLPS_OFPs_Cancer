#
#
import pandas as pd
import numpy as np
from scipy.stats import f_oneway, friedmanchisquare, wilcoxon
from itertools import combinations

LLPS_scaffolds = pd.read_csv("LLPS_scaffolds.csv", header=None)[0].to_list()
LLPS_reg = pd.read_csv("LLPS_regulators.csv", header=None)[0].to_list()
LLPS_clients = pd.read_csv("LLPS_clients.csv", header=None)[0].to_list()
census = pd.read_csv("COSMIC_Census.csv", header=None)[0].to_list()

Ca_scaffolds = list(set(LLPS_scaffolds) & set(census))
Ca_reg = list(set(LLPS_reg) & set(census))
Ca_clients = list(set(LLPS_clients) & set(census))
Ca_rest = list(((set(census) - set(LLPS_scaffolds)) - set(LLPS_reg)) - set(LLPS_clients))

#ST1:
t1 = ["GO:0006306", "GO:0080111", "GO:0030592", "GO:0045006"]
t2 = ["GO:0071103", "GO:0034728", "GO:0006338", "GO:0070828", "GO:0070827", "GO:0006333"]
t3 = ["GO:0006281", "GO:0042769"]
t4 = ["GO:0006310"]
t5 = ["GO:0006260", "GO:0007059"]
#ST2:
t6 = ["GO:0006351", "GO:0006355", "GO:0032922", "GO:0040029", "GO:1900095", "GO:0060968"]
t7 = ["GO:0010608", "GO:0006397", "GO:2000235", "GO:2000232", "GO:0048024", "GO:1902796", "GO:0006412", "GO:0048024", "GO:0006402", "GO:0061013"]
t8 = ["GO:0051604", "GO:0006457", "GO:1903317"]
t9 = ["GO:0044257", "GO:1903362", "GO:0006508", "GO:0031647"]
t10= ["GO:1902414", "GO:0071692", "GO:1905719", "GO:0033365", "GO:0034394", "GO:0051668", "GO:0032507", "GO:0045184", "GO:0032880"]
#ST3:
t11 = ["GO:0006468", "GO:0006470", "GO:0006486", "GO:0006517", "GO:0006479", "GO:0006482", "GO:0043543", "GO:0035601", "GO:0070647", "GO:0006497", "GO:0051697", "GO:0006486", "GO:0006517", "GO:0018126", "GO:0043687", "GO:0032020", "GO:0006471", "GO:0016570"]
t12 = ["GO:0043393", "GO:0051101", "GO:1900130", "GO:1905214", "GO:0035561", "GO:0033341", "GO:1905853", "GO:1905696"]
t13 = ["GO:0065003", "GO:0032984", "GO:0043254", "GO:0043244", "GO:0071823"]
t14 = ["GO:0050790"]
#ST4:
t15 = ["GO:0007166", "GO:0007186"]
t16 = ["GO:0035556"]
t17 = ["GO:0030036", "GO:0045104", "GO:0000226", "GO:0030865", "GO:0032185", "GO:1990933", "GO:0030952"]
t18 = ["GO:0098609", "GO:0033627", "GO:0007160", "GO:0034446", "GO:0006931", "GO:0022407", "GO:0033628", "GO:0001952", "GO:1900024", "GO:1904235"]
#ST5:
t19 = ["GO:0098739", "GO:0140115", "GO:0006897", "GO:0006909", "GO:0006887"]
t20 = ["GO:0046907", "GO:0032386", "GO:0098927", "GO:0048193", "GO:0099075"]
t21 = ["GO:0017144", "GO:0043171", "GO:0006749", "GO:0002936", "GO:0005975", "GO:0006629", "GO:0042445", "GO:0072593", "GO:0044281", "GO:0019674", "GO:0006739", "GO:0006734", "GO:0046034"]

GO_slim = pd.read_csv("../GO_of_non-LLPS_proteins/GOslim-annotations-HumanProteome.tsv", sep="\t")

# -------------------------------
# ANALYSIS OF THE LLPS SETS
# -------------------------------

t_overlap = {}
group_names = ["LLPS_scaffolds", "LLPS_regulators", "LLPS_clients", "COSMIC_Census"]
i = 0
for group in [LLPS_scaffolds, LLPS_reg, LLPS_clients, census]:
    overlaps = []
    group_GO_slim = GO_slim[GO_slim["GENE PRODUCT ID"].isin(group)]
    for t in [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21]:
        t_proteins = group_GO_slim[group_GO_slim["GO TERM"].isin(t)]
        uniprots = t_proteins["GENE PRODUCT ID"].drop_duplicates().to_list()
        overlaps.append(len(uniprots)/len(group))
    t_overlap[ group_names[i] ] = overlaps
    i += 1

toolkits = pd.DataFrame().from_dict(t_overlap)
toolkits.index = ["1. DNA chem. mod.", "2. DNA struct. org.", "3. DNA damage repair", "4. DNA recomb.", \
                  "5. DNA repl.", "6. Transcr, gene expr. reg.", "7. mRNA proc. transl. degr.", \
                  "8. Prot. matur, folding", "9. Altering prot. stab.", "10. Altering-maint. prot. loc.", \
                  "11. Prot. PTM", "12. Modul. macromol. intera.", "13. Mol. (dis)assembly", \
                  "14. Reg. cat. act.", "15. Cell surf. rec. signal.", "16. Intracell. signal. transd.", \
                  "17. Cytoskel. org.", "18. Cell adhesion", "19. Transp. across plasmamembr.", \
                  "20. Transp. inside the cell", "21. Metabolism"]

toolkits["LLPS_regulators"]["1. DNA chem. mod."] += 0.000001
print(toolkits)

# -------------------------------
# Friedman test for repeated samples
# -------------------------------

stat, p = friedmanchisquare(toolkits["LLPS_scaffolds"].to_list(), \
                            toolkits["LLPS_regulators"].to_list(), \
                            toolkits["LLPS_clients"].to_list(), \
                            toolkits["COSMIC_Census"].to_list()  )
print()
print("stat =", stat)
print("p =", p)
print()

# -------------------------------
# Pairwise Wilcoxon tests with Bonferroni correction
# -------------------------------

groups = toolkits.columns.tolist()
results = []

# All pairwise combinations of group columns
for g1, g2 in combinations(groups, 2):
    stat, p = wilcoxon(toolkits[g1], toolkits[g2], zero_method="zsplit", method='exact')
    results.append({
        'Group 1': g1,
        'Group 2': g2,
        'Statistic': stat,
        'Raw p-value': p
    })

# Bonferroni correction
num_tests = len(results)
for res in results:
    corrected_p = min(res['Raw p-value'] * num_tests, 1.0)
    res['Bonferroni corrected p'] = corrected_p
    res['Significant'] = corrected_p < 0.05

# Convert to DataFrame and display
posthoc_df = pd.DataFrame(results)
print(posthoc_df)
print("\n")



# -------------------------------
# ANALYSIS OF THE CANCER-ASSOCIATED SETS
# -------------------------------


t_overlap = {}
group_names = ["Cancer-assoc. LLPS scaffolds", "Cancer-assoc. LLPS regulators", "Cancer-assoc. LLPS clients", "Rest of COSMIC Census"]
i = 0
for group in [Ca_scaffolds, Ca_reg, Ca_clients, Ca_rest]:
    overlaps = []
    group_GO_slim = GO_slim[GO_slim["GENE PRODUCT ID"].isin(group)]
    for t in [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21]:
        t_proteins = group_GO_slim[group_GO_slim["GO TERM"].isin(t)]
        uniprots = t_proteins["GENE PRODUCT ID"].drop_duplicates().to_list()
        overlaps.append(len(uniprots)/len(group))
    t_overlap[ group_names[i] ] = overlaps
    i += 1

toolkits = pd.DataFrame().from_dict(t_overlap)
toolkits.index = ["1. DNA chem. mod.", "2. DNA struct. org.", "3. DNA damage repair", "4. DNA recomb.", \
                  "5. DNA repl.", "6. Transcr, gene expr. reg.", "7. mRNA proc. transl. degr.", \
                  "8. Prot. matur, folding", "9. Altering prot. stab.", "10. Altering-maint. prot. loc.", \
                  "11. Prot. PTM", "12. Modul. macromol. intera.", "13. Mol. (dis)assembly", \
                  "14. Reg. cat. act.", "15. Cell surf. rec. signal.", "16. Intracell. signal. transd.", \
                  "17. Cytoskel. org.", "18. Cell adhesion", "19. Transp. across plasmamembr.", \
                  "20. Transp. inside the cell", "21. Metabolism"]

toolkits["Cancer-assoc. LLPS regulators"]["1. DNA chem. mod."] += 0.000001
toolkits["Cancer-assoc. LLPS clients"]["1. DNA chem. mod."] += 0.000002
print(toolkits)

# -------------------------------
# Friedman test for repeated samples
# -------------------------------

stat, p = friedmanchisquare(toolkits["Cancer-assoc. LLPS scaffolds"].to_list(), \
                            toolkits["Cancer-assoc. LLPS regulators"].to_list(), \
                            toolkits["Cancer-assoc. LLPS clients"].to_list(), \
                            toolkits["Rest of COSMIC Census"].to_list() )
print()
print("stat =", stat)
print("p =", p)
print()

# -------------------------------
# Pairwise Wilcoxon tests with Bonferroni correction
# -------------------------------

groups = toolkits.columns.tolist()
results = []

# All pairwise combinations of group columns
for g1, g2 in combinations(groups, 2):
    stat, p = wilcoxon(toolkits[g1], toolkits[g2], zero_method='zsplit', method='exact')
    results.append({
        'Group 1': g1,
        'Group 2': g2,
        'Statistic': stat,
        'Raw p-value': p
    })

# Bonferroni correction
num_tests = len(results)
for res in results:
    corrected_p = min(res['Raw p-value'] * num_tests, 1.0)
    res['Bonferroni corrected p'] = corrected_p
    res['Significant'] = corrected_p < 0.05

# Convert to DataFrame and display
posthoc_df = pd.DataFrame(results)
print(posthoc_df)
