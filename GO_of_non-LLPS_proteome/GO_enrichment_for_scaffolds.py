#
#
import pandas as pd

df_5 = pd.read_csv("UniProt_annotation_score_5_either_human_proteome_or_reviewed.csv", header=None)
annotation_5 = df_5[0].to_list()
df_go = pd.read_csv("GOslim-annotations-HumanProteome.tsv", sep='\t')
df_go_unip = df_go[df_go["GENE PRODUCT DB"] == "UniProtKB"]
df_query = df_go_unip[["GENE PRODUCT ID", "GO TERM"]]
df_query_5 = df_query[df_query["GENE PRODUCT ID"].isin(annotation_5) ]
df_GOs = df_query.drop_duplicates()

bg_counts = df_GOs["GO TERM"].value_counts().to_dict()
print(df_GOs["GO TERM"].value_counts() )

df_scaffolds = pd.read_csv("LLPS_scaffolds.csv", header=None)
scaffolds = df_scaffolds[0].to_list()

df_scaffolds2 = pd.read_csv("PhaSepDB2_human_entries.csv", header=None)
scaffolds += df_scaffolds2[0].to_list()

df_GOs_scaffolds = df_GOs[df_GOs["GENE PRODUCT ID"].isin(scaffolds) ]

scaffold_counts = df_GOs_scaffolds["GO TERM"].value_counts().to_dict()
print(df_GOs_scaffolds["GO TERM"].value_counts() )
scaffold_bg_ratio = len(scaffolds)/len(annotation_5)
print()
print("GO SCAFFOLDS BACKGROUND EXPECTED FOLD_ENRICHMENT")
count_greater_prot = 0
for go in scaffold_counts.keys():
    print(go, scaffold_counts[go], bg_counts[go], scaffold_bg_ratio*bg_counts[go],  scaffold_counts[go]/(scaffold_bg_ratio*bg_counts[go]) )
