import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import os
import torch
import math
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#########################################
data = pd.read_excel("Hallmarks_enrichments.xlsx", sheet_name="without_tissues2", header=0)
# file=data.iloc[1:10,:]
# file=data.iloc[9:18,:]
data=data.iloc[19:30,:]
xlabels = ["LLPS scaffolds (35)", "Regulators (36)", "Clients (226)"]
characteristics=[]
dataframes_list = []
for n in xlabels:
    new_df = data.filter(["characteristics", str(n) + " P_value(two_sided)",
                          str(n) + " Fold enrichment"])
    new_df["label"] = n
    new_df.columns = ["characteristics", "P_value", "enrichment", "label"]
    ###### Replacing "nan" fold enrichment values with -1 so then we can filter for - values
    new_df["enrichment"] = new_df["enrichment"].fillna(-1)
    ###### Replacing "the NOT APPLICABLE" P_values with 0
    for index, row in new_df.iterrows():
        if row["P_value"] == "Not applicable":
            new_df.loc[index, "-P_value_log10"] = 0
        else:
            new_df.loc[index, "-P_value_log10"] = -(math.log10(row["P_value"]))
    new_df = new_df[::-1]
    dataframes_list.append(new_df)
### Concatenate dataframes and assigning new index
Concatenated_dataframe = pd.concat(dataframes_list, ignore_index=True)
# print(Concatenated_dataframe)
# exit()
size_P_value = Concatenated_dataframe["-P_value_log10"].values
# size_P_value = np.where(size_P_value < 1, 0, size_P_value)
# print(len(size_P_value))
# exit()
Concatenated_dataframe["R"]=size_P_value / size_P_value.max() / 2

# Concatenated_dataframe.set_index("characteristics", inplace = True)
# print(Concatenated_dataframe.loc["splice_site_mutations"])


data=Concatenated_dataframe[(Concatenated_dataframe["characteristics"]=="angiogenesis")|(Concatenated_dataframe["characteristics"]=="replicative immortality")|(Concatenated_dataframe["characteristics"]=="energetics")|
(Concatenated_dataframe["characteristics"]=="immune evasion")|(Concatenated_dataframe["characteristics"]=="apoptosis")|(Concatenated_dataframe["characteristics"]=="genomic instability")|(Concatenated_dataframe["characteristics"]=="metastasis")|(Concatenated_dataframe["characteristics"]=="proliferative signaling")|(Concatenated_dataframe["characteristics"]=="growth suppression")|(Concatenated_dataframe["characteristics"]=="inflammation")]
# data=Concatenated_dataframe[(Concatenated_dataframe["characteristics"]=="somatic mutations")|(Concatenated_dataframe["characteristics"]=="germline mutations")|(Concatenated_dataframe["characteristics"]=="role TSG")|(Concatenated_dataframe["characteristics"]=="role oncogene")|
#                             (Concatenated_dataframe["characteristics"]=="role fusion protein")|(Concatenated_dataframe["characteristics"]=="genetics dominant")|(Concatenated_dataframe["characteristics"]=="genetics recessive")|(Concatenated_dataframe["characteristics"]=="actionability level")|(Concatenated_dataframe["characteristics"]=="actionability level")]
# data=Concatenated_dataframe[(Concatenated_dataframe["characteristics"]=="resistance_mutations")|(Concatenated_dataframe["characteristics"]=="fusion")|(Concatenated_dataframe["characteristics"]=="translocation")|(Concatenated_dataframe["characteristics"]=="missense_mutations")|(Concatenated_dataframe["characteristics"]=="splice_site_mutations")|
#                             (Concatenated_dataframe["characteristics"]=="nonsense_mutations")|(Concatenated_dataframe["characteristics"]=="frameshift_mutations")|(Concatenated_dataframe["characteristics"]=="large deletions")|(Concatenated_dataframe["characteristics"]=="amplification")|(Concatenated_dataframe["characteristics"]=="other_mutations")]

data = data.reset_index(drop=True)

ylabels = data["characteristics"].unique().tolist()
xlabels = data["label"].unique().tolist()
color_enrichment = data["enrichment"].values
xn = len(xlabels)
yn = len(ylabels)

# color_enrichment = color_enrichment.clip(min=-1, max=2)

fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlim(-0.5, xn - 0.5)
ax.set_ylim(-0.5, yn - 0.5)
ax.set(xticks=np.arange(xn), yticks=np.arange(yn))
ax.set_xticklabels(xlabels, rotation=90, fontsize=12)
ax.set_yticklabels(ylabels, rotation='horizontal', fontsize=12)
ax.set_xticks(np.arange(xn) - 0.5, minor=True)
ax.set_yticks(np.arange(yn) - 0.5, minor=True)
ax.grid(which='minor')
ax.set_aspect("equal", "box")


# size_P_value = np.where(size_P_value >0.1, 0, size_P_value)
R=data["R"].values
circle = [plt.Circle((xlabels.index(data.loc[i, "label"]),ylabels.index(data.loc[i, "characteristics"])),radius=r) for i, r in enumerate(R)]
col = PatchCollection(circle, array=color_enrichment, cmap=plt.cm.get_cmap('coolwarm'), clim=(-4, 4))
ax.add_collection(col)
cbar = fig.colorbar(col, label="Fold enrichment", pad=0.05)
ax = cbar.ax
text = ax.yaxis.label
font = matplotlib.font_manager.FontProperties(size=13)
text.set_font_properties(font)

#
class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=orig_handle.width + xdescent,
                             height=orig_handle.height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]

smaxx = 4
smax= 3
smidd = 2
smid=1
sminn = 0
# plotting a test case "two entry legend", tring to get two differently sized, differently colored ellipses
g = [mpatches.Ellipse((), width=smaxx * 5, height=smaxx * 5, color="grey"),
     mpatches.Ellipse((), width=smax * 5, height=smax * 5, color="grey"),
     mpatches.Ellipse((), width=smidd * 5, height=smidd * 5, color="grey"),
     mpatches.Ellipse((), width=smid * 5, height=smid * 5, color="grey"),
     mpatches.Ellipse((), width=sminn * 5, height=sminn * 5, color="grey")]
# Using plt.Circles() ?
g = [plt.Circle((), radius=smaxx * 5, color="grey"),
     plt.Circle((), radius=smax * 5, color="grey"),
     plt.Circle((), radius=smidd * 5, color="grey"),
     plt.Circle((), radius=smid * 5, color="grey"),
     plt.Circle((), radius=sminn * 5, color="grey")]
# Update x-axis labels manually


legend = ax.legend(g, ['≥4', "3",'2', "1", "0"], handler_map={mpatches.Ellipse: HandlerEllipse()}, title="-log10(P_value)",
                   title_fontsize="3", fontsize=12,loc='center', bbox_to_anchor=(5, 0.7, 0.7, -0.5), labelspacing=2.8, frameon=False, markerscale=1)
plt.setp(legend.get_title(), fontsize="large")
title = legend.get_title()
plt.tight_layout(pad=1)
plt.show(block=False)

plt.savefig("2.png")
exit()

