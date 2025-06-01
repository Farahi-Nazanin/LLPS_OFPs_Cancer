import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import plotly.graph_objects as go
import plotly.offline as pyo

import plotly.graph_objs as go



df = pd.read_excel("Hallmarks_enrichments.xlsx", sheet_name="without_tissues2", header=0)
xlabels = ["LLPS scaffolds", "Regulators", "Clients", "PhaSepDB"]
dataframes_list = []
for n in xlabels:
    new_df = df.filter(["characteristics", str(n) + " Fold enrichment"])
    new_df["label"] = n
    new_df.columns = ["characteristics", "enrichment", "label"]
    ###### Replacing "nan" fold enrichment values with -1 so then we can filter for -1 values
    new_df["enrichment"] = new_df["enrichment"].fillna(-1)
    dataframes_list.append(new_df)
### Concatenate dataframes and assigning new index
Concatenated_dataframe = pd.concat(dataframes_list, ignore_index=True)
# print(Concatenated_dataframe)
# exit()
"""""""##############################NEW radar plot""""""
# 
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="polar")
# 
# # theta has 5 different angles, and the first one repeated
# theta = np.arange(len(df) + 1) / float(len(df)) * 2 * np.pi
# # values has the 5 values from 'Col B', with the first element repeated
# values = df['Col B'].values
# values = np.append(values, values[0])
# 
# # draw the polygon and the mark the points for each angle/value combination
# l1, = ax.plot(theta, values, color="C2", marker="o", label="Name of Col B")
# plt.xticks(theta[:-1], df['Col A'], color='grey', size=12)
# ax.tick_params(pad=10) # to increase the distance of the labels to the plot
# # fill the area of the polygon with green and some transparency
# ax.fill(theta, values, 'green', alpha=0.1)
# 
# # plt.legend() # shows the legend, using the label of the line plot (useful when there is more than 1 polygon)
# plt.title("Title")
# plt.show()
# 
"""""




############################# Radar plot with mutation types
# drivers = []
# regulators=[]
# clients=[]
# PhaSepDB=[]
# for index, row in Concatenated_dataframe.iterrows():
#     if row["label"] == "LLPS scaffolds" and ( row["characteristics"] == "fusion/translocation" or row[
#      "characteristics"] == "missense mut" or row["characteristics"] == "splice site mut" or row[
#                 "characteristics"] == "Truncating mut" or row["characteristics"] == "large deletions"
#             or row["characteristics"] == "amplification" or row["characteristics"] == "other mut" ):
#         drivers.append(row["enrichment"])
#     if row["label"] == "Regulators" and (row["characteristics"] == "fusion/translocation" or row[
#      "characteristics"] == "missense mut" or row["characteristics"] == "splice site mut" or row[
#                 "characteristics"] == "Truncating mut" or row["characteristics"] == "large deletions"
#             or row["characteristics"] == "amplification" or row["characteristics"] == "other mut" ):
#         regulators.append(row["enrichment"])
#     if row["label"] == "Clients" and (
#             row["characteristics"] == "fusion/translocation" or row[
#         "characteristics"] == "missense mut" or row["characteristics"] == "splice site mut" or row[
#                 "characteristics"] == "Truncating mut" or row["characteristics"] == "large deletions"
#             or row["characteristics"] == "amplification" or row["characteristics"] == "other mut"):
#         clients.append(row["enrichment"])
#     if row["label"] == "PhaSepDB" and (row["characteristics"] == "fusion/translocation" or row[
#      "characteristics"] == "missense mut" or row["characteristics"] == "splice site mut" or row[
#                 "characteristics"] == "Truncating mut" or row["characteristics"] == "large deletions"
#             or row["characteristics"] == "amplification" or row["characteristics"] == "other mut" ):
#         PhaSepDB.append(row["enrichment"])
#
#
#
#
#
# categories = ["fusion/translocation", "missense mut", "splice site mut","Truncating mut","large deletions","amplification", "other mut"]
# categories = [*categories, categories[0]]
#
# drivers = [*drivers, drivers[0]]
# regulators = [*regulators, regulators[0]]
# clients = [*clients, clients[0]]
# PhaSepDB = [*PhaSepDB, PhaSepDB[0]]
#
#
# fig = go.Figure(
#     data=[
#         go.Scatterpolar(r=drivers, theta=categories, name='LLPS scaffolds', line_width=4),
#         go.Scatterpolar(r=regulators, theta=categories, name='Regulators', line_width=4),
#         go.Scatterpolar(r=clients, theta=categories, name='Clients',line=dict(color="goldenrod"), line_width=4), go.Scatterpolar(r=PhaSepDB, theta=categories, name='PhaSepDB', line_width=4)
#     ],
#     layout=go.Layout(
#         # title=go.layout.Title(text='radar'),
#         polar={'radialaxis': {'visible': True, "tickfont_size": 20}, "angularaxis": {"tickfont_size": 20}},
#         showlegend=True
#     )
# )
# fig.update_polars(radialaxis=dict(range=[-1, 1.5]))
# fig.update_layout(legend={"font":{"size":20}})
# pyo.plot(fig)
#
# exit()
#
#
#
# ############################# 5 category radar plot
drivers = []
regulators=[]
clients=[]
PhaSepDB=[]

for index, row in Concatenated_dataframe.iterrows():
    if row["label"] == "LLPS scaffolds" and (row["characteristics"] == "genetics recessive" or row["characteristics"] == "genetics dominant" or row[
     "characteristics"] == "role TSG" or row["characteristics"] == "role oncogene" or row[
                "characteristics"] == "actionability level"):
        drivers.append(row["enrichment"])
    if row["label"] == "Regulators" and (row["characteristics"] == "genetics recessive" or row["characteristics"] == "genetics dominant" or row[
     "characteristics"] == "role TSG" or row["characteristics"] == "role oncogene" or row[
                "characteristics"] == "actionability level"):
        regulators.append(row["enrichment"])
    if row["label"] == "Clients" and (
            row["characteristics"] == "genetics recessive" or row["characteristics"] == "genetics dominant" or row[
        "characteristics"] == "role TSG" or row["characteristics"] == "role oncogene" or row[
                "characteristics"] == "actionability level"):
        clients.append(row["enrichment"])
    if row["label"] == "PhaSepDB" and (
            row["characteristics"] == "genetics recessive" or row["characteristics"] == "genetics dominant" or row[
        "characteristics"] == "role TSG" or row["characteristics"] == "role oncogene" or row[
                "characteristics"] == "actionability level"):
        PhaSepDB.append(row["enrichment"])
# print(drivers, regulators, clients)
# exit()


## with the order in the excel file
categories = ['Tumor suppressors','Oncogene','Dominant', 'Recessive', 'Actionability']
categories = [*categories, categories[0]]

drivers = [*drivers, drivers[0]]
regulators = [*regulators, regulators[0]]
clients = [*clients, clients[0]]
PhaSepDB = [*PhaSepDB, PhaSepDB[0]]


fig = go.Figure(
    data=[
        go.Scatterpolar(r=drivers, theta=categories, name='LLPS scaffolds', line_width=4, mode = 'lines+markers'),
        go.Scatterpolar(r=regulators, theta=categories, name='Regulators',line_width=4, mode = 'lines+markers'),
        go.Scatterpolar(r=clients, theta=categories, name='Clients', line=dict(color="goldenrod"), line_width=4, mode = 'lines+markers'), go.Scatterpolar(r=PhaSepDB, theta=categories, name='PhaSepDB',line_width=4, mode = 'lines+markers'),
    ],
    layout=go.Layout(
        # title=go.layout.Title(text='radar'),
        polar={'radialaxis': {"range":[-1, 1.5],'visible': True,"tickfont_size":20, "gridcolor":"black"}, "angularaxis":{"tickfont_size":20}},
        showlegend=True, plot_bgcolor = "rgba(0,0,0,0)", paper_bgcolor = "rgba(0,0,0,0)"
    )
)
fig.update_polars(radialaxis=dict(range=[-1, 1.5]))
fig.update_layout(legend={"font":{"size":20}})
pyo.plot(fig)

exit()

