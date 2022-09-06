#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"

##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/6_SAMPLE_FILTRATION"
subprocess.run(cmd, shell=True)


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/5_MICROGLIA_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_Filtered.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()


filtered_rows = []

for row in range(0, len(body)):
    if (float(body[row][header.index("WM_Percentage")]) >= 0.5) or (float(body[row][header.index("NAWM_Percentage")]) >= 0.5) or (float(body[row][header.index("AL_Percentage")]) >= 0.5) or (float(body[row][header.index("CA_Percentage")]) >= 0.5):
        filtered_rows.append(body[row])

fvfile = []
fvfile.append(header)
for row in range(0, len(filtered_rows)):
    fvfile.append(filtered_rows[row])


f = open(f"{OUTPUT_POSTPROCESSING_PATH}/6_SAMPLE_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()

