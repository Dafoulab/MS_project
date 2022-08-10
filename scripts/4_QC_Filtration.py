##############################
##############################
### MS GSM4114296 ANALYSIS ###
##############################
##############################

import csv
import subprocess


OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/postprocessing"

##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/4_QC_FILTRATION"
subprocess.run(cmd, shell=True)


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/3_OC_ANNOTATED/MS_GSM4114296_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_Annotated.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()

filtered_header = []
filtered_body = []

filtered_header.append(header)


rnaedit_index = header.index("RNAEDIT")
homopolymer_index = header.index("Homopolymer")
distance_to_junction_index = header.index("DJ")
entropy_index = header.index("Entropy")
repeat_family_index = header.index("Repeat_Family")
dbsnp_index = header.index("dbSNP_rsID")
pseudogene_hugo_index = header.index("Pseudogene_Hugo")
pseudogene_transcript_index = header.index("Pseudogene_Transcript")
sequence_ontology_index = header.index("Sequence_Ontology")
gnomad3_af_index = header.index("gnomAD3_AF")


pseudogene_name_mappings = ["pseudogene", "processed_pseudogene", "polymorphic_pseudogene", "retrotransposed", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "IG_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "Mt_tRNA_pseudogene", "tRNA_pseudogene", "snoRNA_pseudogene", "snRNA_pseudogene", "scRNA_pseudogene", "rRNA_pseudogene", "misc_RNA_pseudogene", "miRNA_pseudogene"]


for row in range(0, len(body)):
    ### If site present in REDIPortal Database and site not in Pseudogenes then keep them ###
    if (body[row][rnaedit_index] != "NA") and (body[row][sequence_ontology_index] not in pseudogene_name_mappings):
        filtered_body.append(body[row])
    else:
        ### If site in ALU regions and site NOT in dbSNP or NOT in Pseudogenes or gnomAD3_AF <= 0.05, then keep them ###
        if body[row][repeat_family_index] == "Alu":
            if (body[row][dbsnp_index] == "NA") and (body[row][pseudogene_hugo_index] == "NA") and (body[row][pseudogene_transcript_index] == "NA") and (body[row][sequence_ontology_index] not in pseudogene_name_mappings):
               if (body[row][gnomad3_af_index] == "NA") or (float(body[row][gnomad3_af_index]) <= 0.05):
                    filtered_body.append(body[row])
        ### If site not in ALU regions and site not in REDIPortal Database and Pseudogene regions ###
        elif (body[row][repeat_family_index] != "Simple_repeat") and (body[row][repeat_family_index] != "Low_complexity") and (body[row][dbsnp_index] == "NA") and (body[row][pseudogene_hugo_index] == "NA") and (body[row][pseudogene_transcript_index] == "NA") and (body[row][sequence_ontology_index] not in pseudogene_name_mappings) and (int(body[row][homopolymer_index]) == 0) and (int(body[row][distance_to_junction_index]) > 4):
            if (body[row][gnomad3_af_index] == "NA") or (float(body[row][gnomad3_af_index]) <= 0.05):
                filtered_body.append(body[row])



fvfile = []
fvfile.append(filtered_header[0])
for row in range(0, len(filtered_body)):
    fvfile.append(filtered_body[row])


f = open(f"{OUTPUT_POSTPROCESSING_PATH}/4_QC_FILTRATION/MS_GSM4114296_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_Filtered.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()








