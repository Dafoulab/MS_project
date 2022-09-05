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
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/BARPLOT"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/ENRICHED_GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/BARPLOT"
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


########################
### READ WM DE GENES ###
########################

WM_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GENES/MS_GSE138614_WM_PATHWAY_GENES.csv",index_col=0)
WM_DE_GENES = WM_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=WM_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="WM_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

WM_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

WM_KEGG_DATA = WM_KEGG_DATA[WM_KEGG_DATA["P-value"] <= 0.05]

WM_KEGG_DATA = WM_KEGG_DATA.sort_values("P-value", ascending=False)

WM_KEGG_PATHWAY_TERM = WM_KEGG_DATA["Term"].tolist()

WM_KEGG_P_VALUES = WM_KEGG_DATA["P-value"].tolist()

WM_KEGG_ADJUSTED_P_VALUES = WM_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(WM_KEGG_ADJUSTED_P_VALUES)):
    if WM_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        WM_KEGG_PATHWAY_ENRICHED_GENES = WM_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = WM_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/ENRICHED_GENES/WM_{term_name}_GENES.csv", "w")
        for gene in WM_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


WM_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in WM_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in WM_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(WM_KEGG_PATHWAY_TERM, WM_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF WM GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS/BARPLOT/MS_GSE138614_WM_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM")











##########################
### READ NAWM DE GENES ###
##########################

NAWM_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GENES/MS_GSE138614_NAWM_PATHWAY_GENES.csv",index_col=0)
NAWM_DE_GENES = NAWM_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=NAWM_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="NAWM_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

NAWM_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

NAWM_KEGG_DATA = NAWM_KEGG_DATA[NAWM_KEGG_DATA["P-value"] <= 0.05]

NAWM_KEGG_DATA = NAWM_KEGG_DATA.sort_values("P-value", ascending=False)

NAWM_KEGG_PATHWAY_TERM = NAWM_KEGG_DATA["Term"].tolist()

NAWM_KEGG_P_VALUES = NAWM_KEGG_DATA["P-value"].tolist()

NAWM_KEGG_ADJUSTED_P_VALUES = NAWM_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(NAWM_KEGG_ADJUSTED_P_VALUES)):
    if NAWM_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        NAWM_KEGG_PATHWAY_ENRICHED_GENES = NAWM_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = NAWM_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/ENRICHED_GENES/NAWM_{term_name}_GENES.csv", "w")
        for gene in NAWM_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


NAWM_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in NAWM_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in NAWM_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(NAWM_KEGG_PATHWAY_TERM, NAWM_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF NAWM GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS/BARPLOT/MS_GSE138614_NAWM_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM")



















########################
### READ AL DE GENES ###
########################

AL_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GENES/MS_GSE138614_AL_PATHWAY_GENES.csv",index_col=0)
AL_DE_GENES = AL_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=AL_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="AL_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

AL_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

AL_KEGG_DATA = AL_KEGG_DATA[AL_KEGG_DATA["P-value"] <= 0.05]

AL_KEGG_DATA = AL_KEGG_DATA.sort_values("P-value", ascending=False)

AL_KEGG_PATHWAY_TERM = AL_KEGG_DATA["Term"].tolist()

AL_KEGG_P_VALUES = AL_KEGG_DATA["P-value"].tolist()

AL_KEGG_ADJUSTED_P_VALUES = AL_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(AL_KEGG_ADJUSTED_P_VALUES)):
    if AL_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        AL_KEGG_PATHWAY_ENRICHED_GENES = AL_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = AL_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/ENRICHED_GENES/AL_{term_name}_GENES.csv", "w")
        for gene in AL_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


AL_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in AL_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in AL_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(AL_KEGG_PATHWAY_TERM, AL_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF AL GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS/BARPLOT/MS_GSE138614_AL_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL")


















########################
### READ CA DE GENES ###
########################

CA_DE_GENES = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GENES/MS_GSE138614_CA_PATHWAY_GENES.csv",index_col=0)
CA_DE_GENES = CA_DE_GENES.index.tolist()


##############################################################
### RETRIEVE KEGG PATHWAY ANALYSIS RESULTS THROUGH ENRICHR ###
##############################################################

ENRICHR_KEGG = gseapy.enrichr(
    gene_list=CA_DE_GENES,
    gene_sets=["KEGG_2021_Human"],
    organism="Human", 
    description="CA_DE_GENES",
    outdir=f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/ENRICHR_KEGG_DATA",
    cutoff=0.05,
    no_plot=True
)

ENRICHR_KEGG.results.head()

######################
### READ KEGG DATA ###
######################

CA_KEGG_DATA = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/ENRICHR_KEGG_DATA/KEGG_2021_Human.Human.enrichr.reports.txt", sep="\t")

CA_KEGG_DATA = CA_KEGG_DATA[CA_KEGG_DATA["P-value"] <= 0.05]

CA_KEGG_DATA = CA_KEGG_DATA.sort_values("P-value", ascending=False)

CA_KEGG_PATHWAY_TERM = CA_KEGG_DATA["Term"].tolist()

CA_KEGG_P_VALUES = CA_KEGG_DATA["P-value"].tolist()

CA_KEGG_ADJUSTED_P_VALUES = CA_KEGG_DATA["Adjusted P-value"].tolist()

for i in range(0, len(CA_KEGG_ADJUSTED_P_VALUES)):
    if CA_KEGG_ADJUSTED_P_VALUES[i] <= 0.05:
        CA_KEGG_PATHWAY_ENRICHED_GENES = CA_KEGG_DATA.iloc[i]["Genes"].split(";")
        term_name = CA_KEGG_DATA.iloc[i]["Term"]
        f = open(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/ENRICHED_GENES/CA_{term_name}_GENES.csv", "w")
        for gene in CA_KEGG_PATHWAY_ENRICHED_GENES:
            print(gene, file=f)
        f.close()


CA_KEGG_P_VALUES_LOG_10_TRANSFORMED = [-log(y,10) for y in CA_KEGG_P_VALUES]

bar_color = [("darkblue" if adj_p_value <= 0.05 else "darkviolet") for adj_p_value in CA_KEGG_ADJUSTED_P_VALUES]

plt.style.use('ggplot')

plt.barh(CA_KEGG_PATHWAY_TERM, CA_KEGG_P_VALUES_LOG_10_TRANSFORMED, color=bar_color)

plt.xlabel("-log$_{10}$(P-value)", fontsize=12)

plt.title("KEGG PATHWAY ANALYSIS OF CA GENES", wrap=True)

COLORS = ["darkblue", "darkviolet"]
LABELS = [r"FDR$\leq$0.05", r"FDR$>$0.05"]

# map names to colors
cmap = dict(zip(LABELS, COLORS))

# create the rectangles for the legend
from matplotlib.patches import Patch
patches = [Patch(color=v, label=k) for k, v in cmap.items()]

# add the legend
plt.legend(labels=LABELS, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=12, frameon=False)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS/BARPLOT/MS_GSE138614_CA_KEGG_PATHWAYS.pdf", bbox_inches = 'tight')

plt.show()

delete_empty_folders(f"{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA")


