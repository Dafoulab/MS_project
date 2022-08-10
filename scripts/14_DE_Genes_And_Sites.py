##############################
##############################
### MS GSM4114296 ANALYSIS ###
##############################
##############################

import csv
import subprocess
import pandas as pd

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/postprocessing"


##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/SITES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/SITES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/SITES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GRAPHS"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GENES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/SITES"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GRAPHS"
subprocess.run(cmd, shell=True)




###################################################################
###################################################################
### DE SITES AND GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
###################################################################
###################################################################


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSM4114296_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

# Get the header line
header = vfile[0]

# Get the body
body = vfile[1:]

# Close open connection to file
variant_file.close()

GENE_COLUMN_INDEX = header.index("Gene")
NCBI_GENE_DESCRIPTION_INDEX = header.index("NCBI_Gene_Description")

WMvsAL_T_TEST_INDEX = header.index("WMvsAL_T_TEST_PVALUE")
WMvsCA_T_TEST_INDEX = header.index("WMvsCA_T_TEST_PVALUE")
NAWMvsAL_T_TEST_INDEX = header.index("NAWMvsAL_T_TEST_PVALUE")
NAWMvsCA_T_TEST_INDEX = header.index("NAWMvsCA_T_TEST_PVALUE")

WM_PERCENTAGE_INDEX = header.index("WM_Percentage")
NAWM_PERCENTAGE_INDEX = header.index("NAWM_Percentage")
AL_PERCENTAGE_INDEX = header.index("AL_Percentage")
CA_PERCENTAGE_INDEX = header.index("CA_Percentage")



##########################################################################
### WMvsAL DE SITES AND GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
##########################################################################

WMvsAL_DE_SITES = []

for row in range(0, len(body)):
    if float(body[row][WMvsAL_T_TEST_INDEX]) <= 0.05:
        WMvsAL_DE_SITES.append(body[row])

fvfile = []
fvfile.append(header)
for row in range(0, len(WMvsAL_DE_SITES)):
    fvfile.append(WMvsAL_DE_SITES[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/SITES/MS_GSM4114296_WMvsAL_DE_SITES.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()


WMvsAL_DE_GENES = {}

for row in range(0, len(body)):
    if (body[row][GENE_COLUMN_INDEX] != '' or body[row][GENE_COLUMN_INDEX] != 'NA') and (float(body[row][WMvsAL_T_TEST_INDEX]) <= 0.05) and body[row][GENE_COLUMN_INDEX] not in WMvsAL_DE_GENES:
        WMvsAL_DE_GENES[body[row][GENE_COLUMN_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


WMvsAL_df = pd.DataFrame(pd.Series(WMvsAL_DE_GENES).to_frame("NCBI_Gene_Description"))
WMvsAL_df.reset_index(inplace=True)
WMvsAL_df = WMvsAL_df.rename(columns = {'index':'Gene'})
WMvsAL_df.to_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsAL/GENES/MS_GSM4114296_WMvsAL_DE_GENES.csv", na_rep="none", header=True, index=False)





##########################################################################
### WMvsCA DE SITES AND GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
##########################################################################

WMvsCA_DE_SITES = []

for row in range(0, len(body)):
    if float(body[row][WMvsCA_T_TEST_INDEX]) <= 0.05:
        WMvsCA_DE_SITES.append(body[row])

fvfile = []
fvfile.append(header)
for row in range(0, len(WMvsCA_DE_SITES)):
    fvfile.append(WMvsCA_DE_SITES[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/SITES/MS_GSM4114296_WMvsCA_DE_SITES.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()


WMvsCA_DE_GENES = {}

for row in range(0, len(body)):
    if (body[row][GENE_COLUMN_INDEX] != '' or body[row][GENE_COLUMN_INDEX] != 'NA') and (float(body[row][WMvsCA_T_TEST_INDEX]) <= 0.05) and body[row][GENE_COLUMN_INDEX] not in WMvsCA_DE_GENES:
        WMvsCA_DE_GENES[body[row][GENE_COLUMN_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


WMvsCA_df = pd.DataFrame(pd.Series(WMvsCA_DE_GENES).to_frame("NCBI_Gene_Description"))
WMvsCA_df.reset_index(inplace=True)
WMvsCA_df = WMvsCA_df.rename(columns = {'index':'Gene'})
WMvsCA_df.to_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/WMvsCA/GENES/MS_GSM4114296_WMvsCA_DE_GENES.csv", na_rep="none", header=True, index=False)







############################################################################
### NAWMvsAL DE SITES AND GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
############################################################################

NAWMvsAL_DE_SITES = []

for row in range(0, len(body)):
    if float(body[row][NAWMvsAL_T_TEST_INDEX]) <= 0.05:
        NAWMvsAL_DE_SITES.append(body[row])

fvfile = []
fvfile.append(header)
for row in range(0, len(NAWMvsAL_DE_SITES)):
    fvfile.append(NAWMvsAL_DE_SITES[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/SITES/MS_GSM4114296_NAWMvsAL_DE_SITES.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()


NAWMvsAL_DE_GENES = {}

for row in range(0, len(body)):
    if (body[row][GENE_COLUMN_INDEX] != '' or body[row][GENE_COLUMN_INDEX] != 'NA') and (float(body[row][NAWMvsAL_T_TEST_INDEX]) <= 0.05) and body[row][GENE_COLUMN_INDEX] not in NAWMvsAL_DE_GENES:
        NAWMvsAL_DE_GENES[body[row][GENE_COLUMN_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


NAWMvsAL_df = pd.DataFrame(pd.Series(NAWMvsAL_DE_GENES).to_frame("NCBI_Gene_Description"))
NAWMvsAL_df.reset_index(inplace=True)
NAWMvsAL_df = NAWMvsAL_df.rename(columns = {'index':'Gene'})
NAWMvsAL_df.to_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsAL/GENES/MS_GSM4114296_NAWMvsAL_DE_GENES.csv", na_rep="none", header=True, index=False)




############################################################################
### NAWMvsCA DE SITES AND GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION ###
############################################################################

NAWMvsCA_DE_SITES = []

for row in range(0, len(body)):
    if float(body[row][NAWMvsCA_T_TEST_INDEX]) <= 0.05:
        NAWMvsCA_DE_SITES.append(body[row])

fvfile = []
fvfile.append(header)
for row in range(0, len(NAWMvsCA_DE_SITES)):
    fvfile.append(NAWMvsCA_DE_SITES[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/SITES/MS_GSM4114296_NAWMvsCA_DE_SITES.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()


NAWMvsCA_DE_GENES = {}

for row in range(0, len(body)):
    if (body[row][GENE_COLUMN_INDEX] != '' or body[row][GENE_COLUMN_INDEX] != 'NA') and (float(body[row][NAWMvsCA_T_TEST_INDEX]) <= 0.05) and body[row][GENE_COLUMN_INDEX] not in NAWMvsCA_DE_GENES:
        NAWMvsCA_DE_GENES[body[row][GENE_COLUMN_INDEX]] = body[row][NCBI_GENE_DESCRIPTION_INDEX]


NAWMvsCA_df = pd.DataFrame(pd.Series(NAWMvsCA_DE_GENES).to_frame("NCBI_Gene_Description"))
NAWMvsCA_df.reset_index(inplace=True)
NAWMvsCA_df = NAWMvsCA_df.rename(columns = {'index':'Gene'})
NAWMvsCA_df.to_csv(f"{OUTPUT_POSTPROCESSING_PATH}/14_DE_GENES_AND_SITES/NAWMvsCA/GENES/MS_GSM4114296_NAWMvsCA_DE_GENES.csv", na_rep="none", header=True, index=False)

















