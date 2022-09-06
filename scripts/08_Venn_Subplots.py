#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess
from matplotlib_venn import venn2, venn2_circles
import matplotlib.pyplot as plt

HUMAN_TO_MOUSE_DB_INPUT = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/Human2Mouse/HMD_HumanPhenotype.txt"

#####################################
### Create Human to Mouse Gene DB ###
#####################################


import csv

variant_file = open(HUMAN_TO_MOUSE_DB_INPUT)
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

db_body = vfile[0:]

variant_file.close()


human_to_mouse_genes = {}

for row in range(0, len(db_body)):
    if db_body[row][0] not in human_to_mouse_genes:
        human_to_mouse_genes[db_body[row][0]] = []
        human_to_mouse_genes[db_body[row][0]].append(db_body[row][2])
    else:
        human_to_mouse_genes[db_body[row][0]].append(db_body[row][2])



###############################################################
###############################################################
### VENN DIAGRAMS AFTER QC, MICROGLIA AND SAMPLE FILTRATION ###
###############################################################
###############################################################

for (ms,eae) in [("AL", "CLINICAL"), ("CA", "CLINICAL"), ("WM", "NAIVE"), ("NAWM", "PRECLINICAL")]:
    
    MS_SAMPLE = ms
    EAE_SAMPLE = eae

    OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"
    
    EAE_FILE_INPUT = f"/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/EAE_FILES/EAE_ADAR_{EAE_SAMPLE}.tsv"
    
    WMvsNAWMvsALvsCA_FILE_AFTER_QC_MICROGLIA_SAMPLE_PLUS_T_TEST_INPUT = f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv"

    EAE_vs_MS_COMMON_GENES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/COMMON/GENES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_COMMON_GENES.txt"
    EAE_vs_MS_COMMON_SITES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/COMMON/SITES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_SITES_FROM_COMMON_GENES.tsv"
    
    MS_ONLY_GENES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{MS_SAMPLE}_ONLY/GENES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_{MS_SAMPLE}_ONLY_GENES.txt"
    MS_ONLY_SITES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{MS_SAMPLE}_ONLY/SITES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_{MS_SAMPLE}_ONLY_SITES.tsv"

    EAE_ONLY_GENES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{EAE_SAMPLE}_ONLY/GENES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_{EAE_SAMPLE}_ONLY_GENES.txt"
    EAE_ONLY_SITES_OUTPUT = f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{EAE_SAMPLE}_ONLY/SITES/MS_GSE138614_{MS_SAMPLE}vs{EAE_SAMPLE}_{EAE_SAMPLE}_ONLY_SITES.tsv"

    

    ##################################
    ### Create directory structure ###
    ##################################
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{MS_SAMPLE}_ONLY/GENES"
    subprocess.run(cmd, shell=True)
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{MS_SAMPLE}_ONLY/SITES"
    subprocess.run(cmd, shell=True)
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{EAE_SAMPLE}_ONLY/GENES"
    subprocess.run(cmd, shell=True)
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/{EAE_SAMPLE}_ONLY/SITES"
    subprocess.run(cmd, shell=True)
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/COMMON/GENES"
    subprocess.run(cmd, shell=True)
    cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/{MS_SAMPLE}vs{EAE_SAMPLE}/COMMON/SITES"
    subprocess.run(cmd, shell=True)


    ###############################################
    ##### READ EAE FILE AND FIND UNIQUE GENES #####
    ###############################################

    variant_file = open(EAE_FILE_INPUT)
    read_variant_file = csv.reader(variant_file, delimiter="\t")

    vfile = []

    for row in read_variant_file:
        vfile.append(row)

    EAE_header = vfile[0]
    EAE_body = vfile[1:]

    EAE_GENE_COLUMN_INDEX = EAE_header.index("Gene.refGene")

    EAE_genes = []

    for row in range(0, len(EAE_body)):
        if EAE_body[row][EAE_GENE_COLUMN_INDEX] != '' and EAE_body[row][EAE_GENE_COLUMN_INDEX] not in EAE_genes:
            EAE_genes.append(EAE_body[row][EAE_GENE_COLUMN_INDEX])

    EAE_genes_count = len(EAE_genes)


    ##############################################
    ##### READ MS FILE AND FIND UNIQUE GENES #####
    ##############################################

    variant_file = open(WMvsNAWMvsALvsCA_FILE_AFTER_QC_MICROGLIA_SAMPLE_PLUS_T_TEST_INPUT)
    read_variant_file = csv.reader(variant_file, delimiter="\t")

    vfile = []

    for row in read_variant_file:
        vfile.append(row)

    MS_header = vfile[0]

    MS_body = vfile[1:]

    variant_file.close()

    MS_GENE_COLUMN_INDEX = MS_header.index("Gene")
    MS_PERCENTAGE_INDEX = MS_header.index(f"{ms}_Percentage")

    MS_genes = []

    for row in range(0, len(MS_body)):
        if (MS_body[row][MS_GENE_COLUMN_INDEX] != '' or MS_body[row][MS_GENE_COLUMN_INDEX] != 'NA') and (float(MS_body[row][MS_PERCENTAGE_INDEX]) >= 0.5) and MS_body[row][MS_GENE_COLUMN_INDEX] not in MS_genes:
            MS_genes.append(MS_body[row][MS_GENE_COLUMN_INDEX])

    MS_genes_count = len(MS_genes)


    #####################################
    ##### FIND MSvsEAE COMMON GENES #####
    #####################################

    MS_vs_EAE_common_genes = []

    for gene in MS_genes:
        if gene in human_to_mouse_genes:
            for mouse_gene in human_to_mouse_genes[gene]:
                if (mouse_gene in EAE_genes) and (gene not in MS_vs_EAE_common_genes):
                    MS_vs_EAE_common_genes.append(gene)
                    break

    MS_vs_EAE_common_genes_count = len(MS_vs_EAE_common_genes)

    f = open(EAE_vs_MS_COMMON_GENES_OUTPUT, "w")

    for line in MS_vs_EAE_common_genes:
        print(line, file=f)

    f.close()

    ################################################
    ##### FIND MSvsEAE SITES FROM COMMON GENES #####
    ################################################

    MS_vs_EAE_sites_form_common_genes = []

    for row in range(0, len(MS_body)):
        if MS_body[row][MS_GENE_COLUMN_INDEX] in MS_vs_EAE_common_genes and (float(MS_body[row][MS_PERCENTAGE_INDEX]) >= 0.5):
            MS_vs_EAE_sites_form_common_genes.append(MS_body[row])

    MS_vs_EAE_fvfile = []
    MS_vs_EAE_fvfile.append(MS_header)
    for row in range(0, len(MS_vs_EAE_sites_form_common_genes)):
        MS_vs_EAE_fvfile.append(MS_vs_EAE_sites_form_common_genes[row])

    f = open(EAE_vs_MS_COMMON_SITES_OUTPUT, "w")

    for line in MS_vs_EAE_fvfile:
        print("\t".join(line), file=f)

    f.close()



    ##############################
    ##### FIND MS ONLY GENES #####
    ##############################

    MS_only_genes = []

    for row in range(0, len(MS_body)):
        if (MS_body[row][MS_GENE_COLUMN_INDEX] not in MS_vs_EAE_common_genes) and (MS_body[row][MS_GENE_COLUMN_INDEX] not in MS_only_genes) and (float(MS_body[row][MS_PERCENTAGE_INDEX]) >= 0.5):
            MS_only_genes.append(MS_body[row][MS_GENE_COLUMN_INDEX])

    MS_only_genes_count = len(MS_only_genes)

    f = open(MS_ONLY_GENES_OUTPUT, "w")
    for line in MS_only_genes:
        print(line, file=f)
    f.close()


    #################################################
    ##### FIND MS ONLY SITES FROM MS ONLY GENES #####
    #################################################

    MS_only_lines = []

    for row in range(0, len(MS_body)):
        if MS_body[row][MS_GENE_COLUMN_INDEX] not in MS_vs_EAE_common_genes and (float(MS_body[row][MS_PERCENTAGE_INDEX]) >= 0.5):
            MS_only_lines.append(MS_body[row])

    ms_fvfile = []
    ms_fvfile.append(MS_header)
    for row in range(0, len(MS_only_lines)):
        ms_fvfile.append(MS_only_lines[row])

    f = open(MS_ONLY_SITES_OUTPUT, "w")

    for line in ms_fvfile:
        print("\t".join(line), file=f)

    f.close()


    ###############################
    ##### FIND EAE ONLY GENES #####
    ###############################

    MS_vs_EAE_common_mouse_genes = []

    for gene in MS_genes:
        if gene in human_to_mouse_genes:
            for mouse_gene in human_to_mouse_genes[gene]:
                if (mouse_gene in EAE_genes) and (mouse_gene not in MS_vs_EAE_common_mouse_genes):
                    MS_vs_EAE_common_mouse_genes.append(mouse_gene)
                    break

    MS_vs_EAE_common_mouse_genes_count = len(MS_vs_EAE_common_mouse_genes)


    EAE_only_genes = []

    for gene in EAE_genes:
        if (gene not in MS_vs_EAE_common_mouse_genes) and (gene not in EAE_only_genes):
            EAE_only_genes.append(gene)

    EAE_only_genes_count = len(EAE_only_genes)


    f = open(EAE_ONLY_GENES_OUTPUT, "w")
    for line in EAE_only_genes:
        print(line, file=f)
    f.close()


    ###################################################
    ##### FIND EAE ONLY SITES FROM EAE ONLY SITES #####
    ###################################################

    EAE_only_lines = []

    for row in range(0, len(EAE_body)):
        if EAE_body[row][EAE_GENE_COLUMN_INDEX] not in MS_vs_EAE_common_mouse_genes:
            EAE_only_lines.append(EAE_body[row])

    EAE_fvfile = []
    EAE_fvfile.append(EAE_header)
    for row in range(0, len(EAE_only_lines)):
        EAE_fvfile.append(EAE_only_lines[row])

    f = open(EAE_ONLY_SITES_OUTPUT, "w")

    for line in EAE_fvfile:
        print("\t".join(line), file=f)

    f.close()



    if MS_SAMPLE == "WM" and EAE_SAMPLE == "NAIVE":
        WM_only_genes_count = MS_only_genes_count
        NAIVE_only_genes_count = EAE_only_genes_count
        common_WM_vs_NAIVE_genes_count = MS_vs_EAE_common_genes_count
    elif MS_SAMPLE == "NAWM" and EAE_SAMPLE == "PRECLINICAL":
        NAWM_only_genes_count = MS_only_genes_count
        PRECLINICAL_only_genes_count = EAE_only_genes_count
        common_NAWM_vs_PRECLINICAL_genes_count = MS_vs_EAE_common_genes_count
    elif MS_SAMPLE == "AL" and EAE_SAMPLE == "CLINICAL":
        AL_only_genes_count = MS_only_genes_count
        CLINICAL_only_genes_count = EAE_only_genes_count
        common_AL_vs_CLINICAL_genes_count = MS_vs_EAE_common_genes_count
    elif MS_SAMPLE == "CA" and EAE_SAMPLE == "CLINICAL":
        CA_only_genes_count = MS_only_genes_count
        CLINICAL_only_genes_count = EAE_only_genes_count
        common_CA_vs_CLINICAL_genes_count = MS_vs_EAE_common_genes_count



################################
##### Create Venn Diagrams #####
################################

plt.figure(figsize=(10,10)) 
ax = plt.gca()

colors = ['darkviolet','deepskyblue']
v = venn2(subsets=(WM_only_genes_count, NAIVE_only_genes_count, common_WM_vs_NAIVE_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7)
# set text to defined label id
v.get_label_by_id("10").set_text(f"CONTROL WM\n{WM_only_genes_count}")
v.get_label_by_id("01").set_text(f"EAE NAIVE\n{NAIVE_only_genes_count}")

for text in v.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (WM_only_genes_count, NAIVE_only_genes_count, common_WM_vs_NAIVE_genes_count), linestyle='dashed', linewidth=2, color='k')

plt.title("Number of common and distinct potentially edited genes between CONTROL WM and EAE NAIVE samples", loc='center', wrap=True, fontsize=14)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/WMvsNAIVE/MS_GSE138614_WMvsNAIVE_VENN.pdf")

plt.show()




plt.figure(figsize=(10,10)) 
ax = plt.gca()

colors = ['darkviolet','deepskyblue']
v = venn2(subsets=(NAWM_only_genes_count, PRECLINICAL_only_genes_count, common_NAWM_vs_PRECLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7)
# set text to defined label id
v.get_label_by_id("10").set_text(f"MS NAWM\n{NAWM_only_genes_count}")
v.get_label_by_id("01").set_text(f"EAE PRECLINICAL\n{PRECLINICAL_only_genes_count}")

for text in v.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (NAWM_only_genes_count, PRECLINICAL_only_genes_count, common_NAWM_vs_PRECLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k')

plt.title("Number of common and distinct potentially edited genes between MS NAWM and EAE PRECLINICAL samples", loc='center', wrap=True, fontsize=14)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/NAWMvsPRECLINICAL/MS_GSE138614_NAWMvsPRECLINICAL_VENN.pdf")

plt.show()



plt.figure(figsize=(10,10)) 
ax = plt.gca()

colors = ['darkviolet','deepskyblue']
v = venn2(subsets=(AL_only_genes_count, CLINICAL_only_genes_count, common_AL_vs_CLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7)
# set text to defined label id
v.get_label_by_id("10").set_text(f"MS AL\n{AL_only_genes_count}")
v.get_label_by_id("01").set_text(f"EAE CLINICAL\n{CLINICAL_only_genes_count}")

for text in v.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (AL_only_genes_count, CLINICAL_only_genes_count, common_AL_vs_CLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k')

plt.title("Number of common and distinct potentially edited genes between MS AL and EAE CLINICAL samples", loc='center', wrap=True, fontsize=14)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/ALvsCLINICAL/MS_GSE138614_ALvsCLINICAL_VENN.pdf")

plt.show()





plt.figure(figsize=(10,10)) 
ax = plt.gca()

colors = ['darkviolet','deepskyblue']
v = venn2(subsets=(CA_only_genes_count, CLINICAL_only_genes_count, common_CA_vs_CLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7)
# set text to defined label id
v.get_label_by_id("10").set_text(f"MS CA\n{CA_only_genes_count}")
v.get_label_by_id("01").set_text(f"EAE CLINICAL\n{CLINICAL_only_genes_count}")

for text in v.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (CA_only_genes_count, CLINICAL_only_genes_count, common_CA_vs_CLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k')

plt.title("Number of common and distinct potentially edited genes between MS CA and EAE CLINICAL samples", loc='center', wrap=True, fontsize=14)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/CAvsCLINICAL/MS_GSE138614_CAvsCLINICAL_VENN.pdf")

plt.show()






figure, axes = plt.subplots(2, 2, figsize=(15,15))
colors = ['darkviolet','deepskyblue']


v = venn2(subsets=(WM_only_genes_count, NAIVE_only_genes_count, common_WM_vs_NAIVE_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7, ax=axes[0][0])
# set text to defined label id
v.get_label_by_id("10").set_text(f"CONTROL WM\n{WM_only_genes_count}")
v.get_label_by_id("01").set_text(f"EAE NAIVE\n{NAIVE_only_genes_count}")

for text in v.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (WM_only_genes_count, NAIVE_only_genes_count, common_WM_vs_NAIVE_genes_count), linestyle='dashed', linewidth=2, color='k', ax=axes[0][0])


c = venn2(subsets=(NAWM_only_genes_count, PRECLINICAL_only_genes_count, common_NAWM_vs_PRECLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7, ax=axes[0][1])
# set text to defined label id
c.get_label_by_id("10").set_text(f"MS NAWM\n{NAWM_only_genes_count}")
c.get_label_by_id("01").set_text(f"EAE PRECLINICAL\n{PRECLINICAL_only_genes_count}")

for text in c.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (NAWM_only_genes_count, PRECLINICAL_only_genes_count, common_NAWM_vs_PRECLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k', ax=axes[0][1])



x = venn2(subsets=(AL_only_genes_count, CLINICAL_only_genes_count, common_AL_vs_CLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7, ax=axes[1][0])
# set text to defined label id
x.get_label_by_id("10").set_text(f"MS AL\n{AL_only_genes_count}")
x.get_label_by_id("01").set_text(f"EAE CLINICAL\n{CLINICAL_only_genes_count}")

for text in x.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (AL_only_genes_count, CLINICAL_only_genes_count, common_AL_vs_CLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k', ax=axes[1][0])



z = venn2(subsets=(CA_only_genes_count, CLINICAL_only_genes_count, common_CA_vs_CLINICAL_genes_count), set_labels=('',''), set_colors=colors, alpha=0.7, ax=axes[1][1])
# set text to defined label id
z.get_label_by_id("10").set_text(f"MS CA\n{CA_only_genes_count}")
z.get_label_by_id("01").set_text(f"EAE CLINICAL\n{CLINICAL_only_genes_count}")

for text in z.subset_labels:   
    text.set_color('white')   
    text.set_fontsize(12)   
    text.set_fontweight('bold')  

venn2_circles(subsets = (CA_only_genes_count, CLINICAL_only_genes_count, common_CA_vs_CLINICAL_genes_count), linestyle='dashed', linewidth=2, color='k', ax=axes[1][1])


plt.suptitle("Number of common and distinct potentially edited genes between MS and EAE samples", wrap=True, size=18)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/MS_GSE138614_MSvsEAE_VENN.pdf")


plt.show()






























