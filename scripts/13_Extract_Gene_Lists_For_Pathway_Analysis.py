#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess
import pandas as pd

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"


##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GRAPHS"
subprocess.run(cmd, shell=True)


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()




##########
### WM ###
##########

WM_GENES = {}
NCBI_GENE_DESCRIPTION_INDEX = header.index("NCBI_Gene_Description")
WM_PERCENTAGE_INDEX = header.index("WM_Percentage")
GENE_INDEX = header.index("Gene")

for row in range(0, len(body)):
    if (float(body[row][WM_PERCENTAGE_INDEX]) >= 0.5) and (body[row][GENE_INDEX] not in WM_GENES):
        WM_GENES[body[row][GENE_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]

WM_df = pd.DataFrame(pd.Series(WM_GENES).to_frame("NCBI_Gene_Description"))
WM_df.reset_index(inplace=True)
WM_df = WM_df.rename(columns = {'index':'Gene'})
WM_df.to_csv(f'{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/WM/GENES/MS_GSE138614_WM_PATHWAY_GENES.csv', na_rep="none", header=True, index=False)








############
### NAWM ###
############

NAWM_GENES = {}
NCBI_GENE_DESCRIPTION_INDEX = header.index("NCBI_Gene_Description")
NAWM_PERCENTAGE_INDEX = header.index("NAWM_Percentage")
GENE_INDEX = header.index("Gene")

for row in range(0, len(body)):
    if (float(body[row][NAWM_PERCENTAGE_INDEX]) >= 0.5) and (body[row][GENE_INDEX] not in NAWM_GENES):
        NAWM_GENES[body[row][GENE_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


NAWM_df = pd.DataFrame(pd.Series(NAWM_GENES).to_frame("NCBI_Gene_Description"))
NAWM_df.reset_index(inplace=True)
NAWM_df = NAWM_df.rename(columns = {'index':'Gene'})
NAWM_df.to_csv(f'{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/NAWM/GENES/MS_GSE138614_NAWM_PATHWAY_GENES.csv', na_rep="none", header=True, index=False)




##########
### AL ###
##########

AL_GENES = {}
NCBI_GENE_DESCRIPTION_INDEX = header.index("NCBI_Gene_Description")
AL_PERCENTAGE_INDEX = header.index("AL_Percentage")
GENE_INDEX = header.index("Gene")

for row in range(0, len(body)):
    if (float(body[row][AL_PERCENTAGE_INDEX]) >= 0.5) and (body[row][GENE_INDEX] not in AL_GENES):
        AL_GENES[body[row][GENE_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


AL_df = pd.DataFrame(pd.Series(AL_GENES).to_frame("NCBI_Gene_Description"))
AL_df.reset_index(inplace=True)
AL_df = AL_df.rename(columns = {'index':'Gene'})
AL_df.to_csv(f'{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/AL/GENES/MS_GSE138614_AL_PATHWAY_GENES.csv', na_rep="none", header=True, index=False)





##########
### CA ###
##########


CA_GENES = {}
NCBI_GENE_DESCRIPTION_INDEX = header.index("NCBI_Gene_Description")
CA_PERCENTAGE_INDEX = header.index("CA_Percentage")
GENE_INDEX = header.index("Gene")

for row in range(0, len(body)):
    if (float(body[row][CA_PERCENTAGE_INDEX]) >= 0.5) and (body[row][GENE_INDEX] not in CA_GENES):
        CA_GENES[body[row][GENE_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


CA_df = pd.DataFrame(pd.Series(CA_GENES).to_frame("NCBI_Gene_Description"))
CA_df.reset_index(inplace=True)
CA_df = CA_df.rename(columns = {'index':'Gene'})
CA_df.to_csv(f'{OUTPUT_POSTPROCESSING_PATH}/13_PATHWAY_ANALYSIS/CA/GENES/MS_GSE138614_CA_PATHWAY_GENES.csv', na_rep="none", header=True, index=False)


