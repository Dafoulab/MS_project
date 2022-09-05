#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import os
import subprocess
import pandas as pd
import gseapy
import matplotlib.pyplot as plt
from math import log


OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"

##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)

def delete_empty_folders(root):
    deleted = set()
    for current_dir, subdirs, files in os.walk(root, topdown=False):
        still_has_subdirs = any(subdir for subdir in subdirs if os.path.join(current_dir, subdir) not in deleted)
        if not any(files) and not still_has_subdirs:
            os.rmdir(current_dir)
            deleted.add(current_dir)
    return deleted

#################################
### CHECK AVAILABLE DATABASES ###
#################################

enrichr_available_databases = gseapy.get_library_name()
print(enrichr_available_databases)


############################
### READ WMvsCA DE GENES ###
############################

WMvsCA_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GENES/MS_GSE138614_WMvsCA_DE_GENES.csv",index_col=0)
WMvsCA_DE_GENES = WMvsCA_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=WMvsCA_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="WMvsCA_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

WMvsCA_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

WMvsCA_KEGG_DATA = WMvsCA_KEGG_DATA[WMvsCA_KEGG_DATA["P-value"] <= 0.05]

WMvsCA_KEGG_DATA = WMvsCA_KEGG_DATA.sort_values("P-value", ascending=False)

WMvsCA_KEGG_PATHWAY_TERM = WMvsCA_KEGG_DATA["Term"].tolist()

WMvsCA_KEGG_P_VALUES = WMvsCA_KEGG_DATA["P-value"].tolist()

WMvsCA_KEGG_ADJUSTED_P_VALUES = WMvsCA_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(WMvsCA_KEGG_ADJUSTED_P_VALUES)):
    if WMvsCA_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        WMvsCA_KEGG_PATHWAY_ENRICHED_GENES = WMvsCA_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = WMvsCA_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHED_GENES/WMvsCA_{term_name}_GENES.csv", "w")
        for gene in WMvsCA_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


WMvsCA_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in WMvsCA_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in WMvsCA_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(WMvsCA_KEGG_PATHWAY_TERM, WMvsCA_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF DIFFERENTIALLY EDITED WMvsCA GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/BARPLOT/MS_GSE138614_WMvsAL_DE_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA")








############################
### READ WMvsCA DE GENES ###
############################

WMvsCA_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GENES/MS_GSE138614_WMvsCA_DE_GENES.csv",index_col=0)
WMvsCA_DE_GENES = WMvsCA_DE_GENES.index.tolist()


################################################################
### RETRIEVE MSIGDB PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
################################################################

ENRICHR_MSIGDB = gseapy.enrichr(
    gene_list=WMvsCA_DE_GENES,
    gene_sets=["MSigDB_Hallmark_2020"],
    organism="Human", 
    description="WMvsCA_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHR_MSIGDB_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_MSIGDB.results.head()

########################
### READ MSIGDB DATA ###
########################

WMvsCA_MSIGDB_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHR_MSIGDB_DATA/MSigDB_Hallmark_2020.Human.enrichr.reports.txt", sep="\t")

WMvsCA_MSIGDB_DATA = WMvsCA_MSIGDB_DATA[WMvsCA_MSIGDB_DATA["P-value"] <= 0.05]

WMvsCA_MSIGDB_DATA = WMvsCA_MSIGDB_DATA.sort_values("P-value", ascending=False)

WMvsCA_MSIGDB_PATHWAY_TERM = WMvsCA_MSIGDB_DATA["Term"].tolist()

WMvsCA_MSIGDB_P_VALUES = WMvsCA_MSIGDB_DATA["P-value"].tolist()

WMvsCA_MSIGDB_ADJUSTED_P_VALUES = WMvsCA_MSIGDB_DATA["Adjusted P-value"].tolist()

for i in range(0, len(WMvsCA_MSIGDB_ADJUSTED_P_VALUES)):
    if WMvsCA_MSIGDB_ADJUSTED_P_VALUES[i] <= 0.05:
        WMvsCA_MSIGDB_PATHWAY_ENRICHED_GENES = WMvsCA_MSIGDB_DATA.iloc[i]["Genes"].split(";")
        term_name = WMvsCA_MSIGDB_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/ENRICHED_GENES/WMvsCA_{term_name}_GENES.csv", "w")
        for gene in WMvsCA_MSIGDB_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


WMvsCA_MSIGDB_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in WMvsCA_MSIGDB_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in WMvsCA_MSIGDB_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(WMvsCA_MSIGDB_PATHWAY_TERM, WMvsCA_MSIGDB_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("MSigDB Hallmark 2020 PATHWAY ANALYSIS OF DIFFERENTIALLY EDITED WMvsCA GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS/BARPLOT/MS_GSE138614_WMvsAL_DE_MSigDB_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA")































############################
### READ WMvsAL DE GENES ###
############################

WMvsAL_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GENES/MS_GSE138614_WMvsAL_DE_GENES.csv",index_col=0)
WMvsAL_DE_GENES = WMvsAL_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=WMvsAL_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="WMvsAL_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

WMvsAL_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

WMvsAL_KEGG_DATA = WMvsAL_KEGG_DATA[WMvsAL_KEGG_DATA["P-value"] <= 0.05]

WMvsAL_KEGG_DATA = WMvsAL_KEGG_DATA.sort_values("P-value", ascending=False)

WMvsAL_KEGG_PATHWAY_TERM = WMvsAL_KEGG_DATA["Term"].tolist()

WMvsAL_KEGG_P_VALUES = WMvsAL_KEGG_DATA["P-value"].tolist()

WMvsAL_KEGG_ADJUSTED_P_VALUES = WMvsAL_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(WMvsAL_KEGG_ADJUSTED_P_VALUES)):
    if WMvsAL_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        WMvsAL_KEGG_PATHWAY_ENRICHED_GENES = WMvsAL_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = WMvsAL_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/ENRICHED_GENES/WMvsAL_{term_name}_GENES.csv", "w")
        for gene in WMvsAL_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


WMvsAL_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in WMvsAL_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in WMvsAL_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(WMvsAL_KEGG_PATHWAY_TERM, WMvsAL_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF DIFFERENTIALLY EDITED WMvsAL GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS/BARPLOT/MS_GSE138614_WMvsAL_DE_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()


delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL")


















##############################
### READ NAWMvsCA DE GENES ###
##############################

NAWMvsCA_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GENES/MS_GSE138614_NAWMvsCA_DE_GENES.csv",index_col=0)
NAWMvsCA_DE_GENES = NAWMvsCA_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=NAWMvsCA_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="NAWMvsCA_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

NAWMvsCA_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

NAWMvsCA_KEGG_DATA = NAWMvsCA_KEGG_DATA[NAWMvsCA_KEGG_DATA["P-value"] <= 0.05]

NAWMvsCA_KEGG_DATA = NAWMvsCA_KEGG_DATA.sort_values("P-value", ascending=False)

NAWMvsCA_KEGG_PATHWAY_TERM = NAWMvsCA_KEGG_DATA["Term"].tolist()

NAWMvsCA_KEGG_P_VALUES = NAWMvsCA_KEGG_DATA["P-value"].tolist()

NAWMvsCA_KEGG_ADJUSTED_P_VALUES = NAWMvsCA_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(NAWMvsCA_KEGG_ADJUSTED_P_VALUES)):
    if NAWMvsCA_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        NAWMvsCA_KEGG_PATHWAY_ENRICHED_GENES = NAWMvsCA_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = NAWMvsCA_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/ENRICHED_GENES/NAWMvsCA_{term_name}_GENES.csv", "w")
        for gene in NAWMvsCA_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


NAWMvsCA_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in NAWMvsCA_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in NAWMvsCA_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(NAWMvsCA_KEGG_PATHWAY_TERM, NAWMvsCA_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF DIFFERENTIALLY EDITED NAWMvsCA GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS/BARPLOT/MS_GSE138614_NAWMvsCA_DE_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()


delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA")

















##############################
### READ NAWMvsAL DE GENES ###
##############################

NAWMvsAL_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GENES/MS_GSE138614_NAWMvsAL_DE_GENES.csv",index_col=0)
NAWMvsAL_DE_GENES = NAWMvsAL_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=NAWMvsAL_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="NAWMvsAL_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

NAWMvsAL_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

NAWMvsAL_KEGG_DATA = NAWMvsAL_KEGG_DATA[NAWMvsAL_KEGG_DATA["P-value"] <= 0.05]

NAWMvsAL_KEGG_DATA = NAWMvsAL_KEGG_DATA.sort_values("P-value", ascending=False)

NAWMvsAL_KEGG_PATHWAY_TERM = NAWMvsAL_KEGG_DATA["Term"].tolist()

NAWMvsAL_KEGG_P_VALUES = NAWMvsAL_KEGG_DATA["P-value"].tolist()

NAWMvsAL_KEGG_ADJUSTED_P_VALUES = NAWMvsAL_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(NAWMvsAL_KEGG_ADJUSTED_P_VALUES)):
    if NAWMvsAL_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        NAWMvsAL_KEGG_PATHWAY_ENRICHED_GENES = NAWMvsAL_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = NAWMvsAL_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/ENRICHED_GENES/NAWMvsAL_{term_name}_GENES.csv", "w")
        for gene in NAWMvsAL_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


NAWMvsAL_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in NAWMvsAL_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in NAWMvsAL_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(NAWMvsAL_KEGG_PATHWAY_TERM, NAWMvsAL_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF DIFFERENTIALLY EDITED NAWMvsAL GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS/BARPLOT/MS_GSE138614_NAWMvsAL_DE_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()


delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL")


