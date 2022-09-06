#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"

MICROGLIA_ENRICHED_GENES_FILE = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/Microglia/Microglia_Enriched_Genes.csv"


##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/5_MICROGLIA_FILTRATION"
subprocess.run(cmd, shell=True)



variant_file = open(MICROGLIA_ENRICHED_GENES_FILE)
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

microglia_enriched_genes_vfile_header = vfile[0]
microglia_enriched_genes_vfile_body = vfile[1:]


microglia_gene_index = microglia_enriched_genes_vfile_header.index("hgnc_symbol")
microglia_species_index = microglia_enriched_genes_vfile_header.index("Species")


microglia_enriched_genes = []

for row in range(0, len(microglia_enriched_genes_vfile_body)):
    if (microglia_enriched_genes_vfile_body[row][microglia_gene_index] not in microglia_enriched_genes) and (microglia_enriched_genes_vfile_body[row][microglia_species_index] == "human"):
        microglia_enriched_genes.append(microglia_enriched_genes_vfile_body[row][microglia_gene_index])

variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/4_QC_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_Filtered.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

header = vfile[0]
body = vfile[1:]

variant_file.close()

filtered_rows = []

for row in range(0, len(body)):
  if body[row][5] in microglia_enriched_genes:
    filtered_rows.append(body[row])


fvfile = []
fvfile.append(header)
for row in range(0, len(filtered_rows)):
    fvfile.append(filtered_rows[row])


f = open(f"{OUTPUT_POSTPROCESSING_PATH}/5_MICROGLIA_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_Filtered.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()









