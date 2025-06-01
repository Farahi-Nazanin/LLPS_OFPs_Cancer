import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Load your Excel file
df = pd.read_excel("cosmic_determined_drivers.xlsx")

# Preview column names to identify the exact hallmark columns
print(df.columns[29:39])  # Python is 0-indexed, so columns 28–38 are 27–37
# Sum hallmark columns (binary) to get total per gene
df["hallmark_count"] = df.iloc[:, 29:39].sum(axis=1)

# Create a boolean for LLPS driver status
df["is_LLPS_driver"] = df["driver"].str.lower() == "driver"

# Separate data into two groups
llps_hallmarks = df[df["is_LLPS_driver"]]["hallmark_count"]
non_llps_cosmic_genes = df[~df["is_LLPS_driver"]]["hallmark_count"]

# Mann–Whitney U test (directional: LLPS > non-LLPS)
stat, p_value = mannwhitneyu(llps_hallmarks, non_llps_cosmic_genes, alternative='greater')

# Report result
print(f"Mann–Whitney U test:\nU statistic = {stat:.2f}\np-value = {p_value:.4f}")


# === 6. Bin hallmark counts (group ≥5 together) ===
def bin_counts(series):
    binned = series.copy()
    binned[binned >= 5] = 5
    return binned.value_counts().sort_index()


llps_binned = bin_counts(llps_hallmarks)
non_llps_binned = bin_counts(non_llps_cosmic_genes)


# Ensure all bins are present (0–5)
all_bins = range(6)
llps_binned = llps_binned.reindex(all_bins, fill_value=0)
non_llps_binned = non_llps_binned.reindex(all_bins, fill_value=0)

# === 7. Convert to percentage ===
llps_percent = (llps_binned / llps_binned.sum())* 100
non_llps_percent = (non_llps_binned / non_llps_binned.sum() )* 100

# === 8. Plot bar chart ===
labels = ['0', '1', '2', '3', '4', '≥5']
x = range(len(labels))
width = 0.35

fig, ax = plt.subplots()
ax.bar([i + width / 2 for i in x], llps_percent, width, label='Cancer-associated LLPS scaffolds', color='tomato')
ax.bar([i - width / 2 for i in x], non_llps_percent, width, label='Other COSMIC Census genes (684)', color='royalblue')

# Formatting
ax.set_xlabel('Number of assigned hallmarks', fontsize=12)
ax.set_ylabel('Percentage of genes in each category', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(0, 80)
ax.legend()
plt.show()

