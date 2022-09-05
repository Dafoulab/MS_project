#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"


##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/12_ALU_DISTRIBUTION"
subprocess.run(cmd, shell=True)



variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()



WM_PERCENTAGE_INDEX = header.index("WM_Percentage")
NAWM_PERCENTAGE_INDEX = header.index("NAWM_Percentage")
AL_PERCENTAGE_INDEX = header.index("AL_Percentage")
CA_PERCENTAGE_INDEX = header.index("CA_Percentage")

REPEAT_FAMILY_INDEX = header.index("Repeat_Family")
REPEAT_NAME_INDEX = header.index("Repeat_Name")


WM_ALU = {}

for row in range(0, len(body)):
    if (float(body[row][WM_PERCENTAGE_INDEX]) >= 0.5):
        if (body[row][REPEAT_FAMILY_INDEX] == "Alu"):
            if (body[row][REPEAT_NAME_INDEX] not in WM_ALU):
                WM_ALU[body[row][REPEAT_NAME_INDEX]] = 1
            else:
                WM_ALU[body[row][REPEAT_NAME_INDEX]] += 1




NAWM_ALU = {}

for row in range(0, len(body)):
    if (float(body[row][NAWM_PERCENTAGE_INDEX]) >= 0.5):
        if (body[row][REPEAT_FAMILY_INDEX] == "Alu"):
            if (body[row][REPEAT_NAME_INDEX] not in NAWM_ALU):
                NAWM_ALU[body[row][REPEAT_NAME_INDEX]] = 1
            else:
                NAWM_ALU[body[row][REPEAT_NAME_INDEX]] += 1



AL_ALU = {}

for row in range(0, len(body)):
    if (float(body[row][AL_PERCENTAGE_INDEX]) >= 0.5):
        if (body[row][REPEAT_FAMILY_INDEX] == "Alu"):
            if (body[row][REPEAT_NAME_INDEX] not in AL_ALU):
                AL_ALU[body[row][REPEAT_NAME_INDEX]] = 1
            else:
                AL_ALU[body[row][REPEAT_NAME_INDEX]] += 1



CA_ALU = {}

for row in range(0, len(body)):
    if (float(body[row][CA_PERCENTAGE_INDEX]) >= 0.5):
        if (body[row][REPEAT_FAMILY_INDEX] == "Alu"):
            if (body[row][REPEAT_NAME_INDEX] not in CA_ALU):
                CA_ALU[body[row][REPEAT_NAME_INDEX]] = 1
            else:
                CA_ALU[body[row][REPEAT_NAME_INDEX]] += 1



ALL_ALU = {}

for key, value in WM_ALU.items():
    if key not in ALL_ALU:
        ALL_ALU[key] = []

for key, value in NAWM_ALU.items():
    if key not in ALL_ALU:
        ALL_ALU[key] = []

for key, value in AL_ALU.items():
    if key not in ALL_ALU:
        ALL_ALU[key] = []

for key, value in CA_ALU.items():
    if key not in ALL_ALU:
        ALL_ALU[key] = []       


for key in ALL_ALU.keys():
    if key not in WM_ALU:
        ALL_ALU[key].append(0)
    else:
        ALL_ALU[key].append(WM_ALU[key])
    if key not in NAWM_ALU:
        ALL_ALU[key].append(0)
    else:
        ALL_ALU[key].append(NAWM_ALU[key])
    if key not in AL_ALU:
        ALL_ALU[key].append(0)
    else:
        ALL_ALU[key].append(AL_ALU[key])
    if key not in CA_ALU:
        ALL_ALU[key].append(0)
    else:
        ALL_ALU[key].append(CA_ALU[key])





aludf = pd.DataFrame(ALL_ALU)

# Get TOP 10 Columns
cols = [col for col in aludf.columns if aludf[col].min() > 55]

aludf2 = aludf[cols]

aludf2.insert(loc=0, column='index', value=['WM','NAWM','AL','CA'])

ax = aludf2.plot(kind='bar', stacked=True, figsize=(12,12))

plt.title('The distribution of RNA editing events in Top 10 Alu subfamilies', wrap=True, fontsize=16)

plt.xlabel('Sample Type', fontsize=14)
plt.ylabel('Number of Editing Events in Top 10 Alu subfamilies', fontsize=14)

#rotate x-axis labels
plt.xticks(rotation=0)

ax.set_xticklabels(['CONTROL WM', 'MS NAWM', 'MS AL', 'MS CA'])


plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/12_ALU_DISTRIBUTION/MS_GSE138614_ALU_DISTRIBUTION.pdf")

plt.show()
