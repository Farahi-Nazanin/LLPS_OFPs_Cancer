#
#
import pandas as pd
import numpy as np
import random

df_scaffolds_go = pd.read_excel("GO_enrichment_for_scaffolds.xlsx")
df_scaffolds_go = df_scaffolds_go[df_scaffolds_go["FOLD_ENRICHMENT"] > 1.0]
GOs_by_expectancy  = df_scaffolds_go.sort_values("EXPECTED", ascending=False)["GO"].to_list()

#df_go2unip = pd.read_csv("Bg_annot5_w_LLPSscaffold_GOs.csv")
#all_unips = df_go2unip["GENE PRODUCT ID"].drop_duplicates().to_list()
go2unip = {}
unip2go = {}
with open("Bg_annot5_w_LLPSscaffold_GOs.csv") as infile:
    for line in infile:
        line = line.strip()
        if len(line)>0:
            go, unip = line.split(',')
            if go != "GO TERM": #HEADER
                if go not in go2unip.keys():
                    go2unip[go] = []
                go2unip[go].append(unip)
                if unip not in unip2go.keys():
                    unip2go[unip] = []
                unip2go[unip].append(go)

sel_go2unip = {}
for k in range(1,101):
    for enriched_GO in GOs_by_expectancy:
        GOs_currently = 0
        if enriched_GO in sel_go2unip.keys():
            GOs_currently = len(sel_go2unip[enriched_GO])
        if enriched_GO in go2unip and len(go2unip[enriched_GO])>0:
            uniprots = go2unip[enriched_GO]
            while GOs_currently < int(df_scaffolds_go[df_scaffolds_go["GO"] == enriched_GO]["SCAFFOLDS"]):
                protein = uniprots[ random.randint(0,len(uniprots)-1) ]
                uniprots.remove(protein)
                GOs = unip2go[protein]
                for GO in GOs:
                    if GO not in sel_go2unip.keys():
                        sel_go2unip[GO] = []
                    sel_go2unip[GO].append(protein)
                GOs_currently += 1
            #From new sampling rounds, remove the seleced proteins
            for unip in sel_go2unip[enriched_GO]:
                for go in go2unip.keys():
                    if unip in go2unip[go]:
                        if len(go2unip[go]) == 1:
                            go2unip[go] = []
                        else:
                            go2unip[go].remove(unip)
                if unip in unip2go:
                    del unip2go[unip]
                    
    all_unips = []
    for u in sel_go2unip.values():
        for i in u:
            if i not in all_unips:
                all_unips.append(i)
    sampled_unips = np.random.choice(all_unips, 141)
    outfile = open("bg_sample_"+str(k)+".txt", 'w')
    for id in list(sampled_unips):
        outfile.write(id+"\n")
    outfile.close()
