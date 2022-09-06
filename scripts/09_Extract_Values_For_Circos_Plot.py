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
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/9_CIRCOS"
subprocess.run(cmd, shell=True)


###########################################
### Load WMvsNAWMvsALvsCA filtered file ###
###########################################
variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

WMvsNAWMvsALvsCA_header = vfile[0]

WMvsNAWMvsALvsCA_body = vfile[1:]

variant_file.close()

####################################
### Load EAE NAIVE filtered file ###
####################################
variant_file = open("/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/EAE_FILES/EAE_ADAR_NAIVE.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

EAE_NAIVE_header = vfile[0]

EAE_NAIVE_body = vfile[1:]

variant_file.close()

##########################################
### Load EAE PRECLINICAL filtered file ###
##########################################
variant_file = open("/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/EAE_FILES/EAE_ADAR_PRECLINICAL.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

EAE_PRECLINICAL_header = vfile[0]

EAE_PRECLINICAL_body = vfile[1:]

variant_file.close()

#######################################
### Load EAE CLINICAL filtered file ###
#######################################
variant_file = open("/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/resources/EAE_FILES/EAE_ADAR_CLINICAL.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

EAE_CLINICAL_header = vfile[0]

EAE_CLINICAL_body = vfile[1:]

variant_file.close()


########################################
### Load WMvsNAIVE Common Sites file ###
########################################
variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/WMvsNAIVE/COMMON/SITES/MS_GSE138614_WMvsNAIVE_SITES_FROM_COMMON_GENES.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

WMvsNAIVE_header = vfile[0]

WMvsNAIVE_body = vfile[1:]

variant_file.close()


################################################
### Load NAWMvsPRECLINICAL Common Sites file ###
################################################
variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/NAWMvsPRECLINICAL/COMMON/SITES/MS_GSE138614_NAWMvsPRECLINICAL_SITES_FROM_COMMON_GENES.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

NAWMvsPRECLINICAL_header = vfile[0]

NAWMvsPRECLINICAL_body = vfile[1:]

variant_file.close()



###########################################
### Load ALvsCLINICAL Common Sites file ###
###########################################
variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/ALvsCLINICAL/COMMON/SITES/MS_GSE138614_ALvsCLINICAL_SITES_FROM_COMMON_GENES.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

ALvsCLINICAL_header = vfile[0]

ALvsCLINICAL_body = vfile[1:]

variant_file.close()



###########################################
### Load CAvsCLINICAL Common Sites file ###
###########################################
variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/8_VENN_DIAGRAMS/CAvsCLINICAL/COMMON/SITES/MS_GSE138614_CAvsCLINICAL_SITES_FROM_COMMON_GENES.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
  vfile.append(row)

CAvsCLINICAL_header = vfile[0]

CAvsCLINICAL_body = vfile[1:]

variant_file.close()










#######################################################
##### WM AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#######################################################

WM_PERCENTAGE_INDEX = WMvsNAWMvsALvsCA_header.index("WM_Percentage")

WM_total_edited_sites = 0

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
  if float(WMvsNAWMvsALvsCA_body[row][WM_PERCENTAGE_INDEX]) >= 0.5:
    WM_total_edited_sites += 1

print(WM_total_edited_sites)
# 3258

functional_types = {}

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
    if WMvsNAWMvsALvsCA_body[row][7] != '' and WMvsNAWMvsALvsCA_body[row][7] != 'NA' and (float(WMvsNAWMvsALvsCA_body[row][WM_PERCENTAGE_INDEX]) >= 0.5) and WMvsNAWMvsALvsCA_body[row][7] not in functional_types:
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] = 1
    elif WMvsNAWMvsALvsCA_body[row][7] in functional_types and (float(WMvsNAWMvsALvsCA_body[row][WM_PERCENTAGE_INDEX]) >= 0.5):
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] += 1

print(functional_types)

# intronic 2275
# 3'5'UTR 563
# downstream 64
# upstream 6
# lncRNA 346
# exonic 4




#########################################################
##### NAWM AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#########################################################

NAWM_PERCENTAGE_INDEX = WMvsNAWMvsALvsCA_header.index("NAWM_Percentage")

NAWM_total_edited_sites = 0

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
  if float(WMvsNAWMvsALvsCA_body[row][NAWM_PERCENTAGE_INDEX]) >= 0.5:
    NAWM_total_edited_sites += 1

print(NAWM_total_edited_sites)
# 3214

functional_types = {}

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
    if WMvsNAWMvsALvsCA_body[row][7] != '' and WMvsNAWMvsALvsCA_body[row][7] != 'NA' and (float(WMvsNAWMvsALvsCA_body[row][NAWM_PERCENTAGE_INDEX]) >= 0.5) and WMvsNAWMvsALvsCA_body[row][7] not in functional_types:
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] = 1
    elif WMvsNAWMvsALvsCA_body[row][7] in functional_types and (float(WMvsNAWMvsALvsCA_body[row][NAWM_PERCENTAGE_INDEX]) >= 0.5):
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] += 1

print(functional_types)

# intronic 2152
# 3'5'UTR 644
# downstream 72
# upstream 16
# lncRNA 326
# exonic 4




#######################################################
##### AL AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#######################################################

AL_PERCENTAGE_INDEX = WMvsNAWMvsALvsCA_header.index("AL_Percentage")

AL_total_edited_sites = 0

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
  if float(WMvsNAWMvsALvsCA_body[row][AL_PERCENTAGE_INDEX]) >= 0.5:
    AL_total_edited_sites += 1

print(AL_total_edited_sites)
# 2874

functional_types = {}

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
    if WMvsNAWMvsALvsCA_body[row][7] != '' and WMvsNAWMvsALvsCA_body[row][7] != 'NA' and (float(WMvsNAWMvsALvsCA_body[row][AL_PERCENTAGE_INDEX]) >= 0.5) and WMvsNAWMvsALvsCA_body[row][7] not in functional_types:
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] = 1
    elif WMvsNAWMvsALvsCA_body[row][7] in functional_types and (float(WMvsNAWMvsALvsCA_body[row][AL_PERCENTAGE_INDEX]) >= 0.5):
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] += 1

print(functional_types)

# intronic 1355
# 3'5'UTR 989
# downstream 119
# upstream 10
# lncRNA 394
# exonic 3




#######################################################
##### CA AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#######################################################

CA_PERCENTAGE_INDEX = WMvsNAWMvsALvsCA_header.index("CA_Percentage")

CA_total_edited_sites = 0

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
  if float(WMvsNAWMvsALvsCA_body[row][CA_PERCENTAGE_INDEX]) >= 0.5:
    CA_total_edited_sites += 1

print(CA_total_edited_sites)
# 2051

functional_types = {}

for row in range(0, len(WMvsNAWMvsALvsCA_body)):
    if WMvsNAWMvsALvsCA_body[row][7] != '' and WMvsNAWMvsALvsCA_body[row][7] != 'NA' and (float(WMvsNAWMvsALvsCA_body[row][CA_PERCENTAGE_INDEX]) >= 0.5) and WMvsNAWMvsALvsCA_body[row][7] not in functional_types:
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] = 1
    elif WMvsNAWMvsALvsCA_body[row][7] in functional_types and (float(WMvsNAWMvsALvsCA_body[row][CA_PERCENTAGE_INDEX]) >= 0.5):
        functional_types[WMvsNAWMvsALvsCA_body[row][7]] += 1

print(functional_types)

# intronic 1193
# 3'5'UTR 452
# downstream 78
# upstream 16
# lncRNA 306
# exonic 5




#####################
##### EAE NAIVE #####
#####################

EAE_NAIVE_total_edited_sites = len(EAE_NAIVE_body)

print(EAE_NAIVE_total_edited_sites)
# 775

functional_types = {}

for row in range(0, len(EAE_NAIVE_body)):
    if EAE_NAIVE_body[row][5] != '' and EAE_NAIVE_body[row][5] != 'NA' and EAE_NAIVE_body[row][5] not in functional_types:
        functional_types[EAE_NAIVE_body[row][5]] = 1
    elif EAE_NAIVE_body[row][5] in functional_types:
        functional_types[EAE_NAIVE_body[row][5]] += 1

print(functional_types)

# intronic 408
# 3'5'UTR 185
# downstream 32
# upstream 4
# lncRNA 0
# exonic 146




###########################
##### EAE PRECLINICAL #####
###########################

EAE_PRECLINICAL_total_edited_sites = len(EAE_PRECLINICAL_body)

print(EAE_PRECLINICAL_total_edited_sites)
# 743

functional_types = {}

for row in range(0, len(EAE_PRECLINICAL_body)):
    if EAE_PRECLINICAL_body[row][5] != '' and EAE_PRECLINICAL_body[row][5] != 'NA' and EAE_PRECLINICAL_body[row][5] not in functional_types:
        functional_types[EAE_PRECLINICAL_body[row][5]] = 1
    elif EAE_PRECLINICAL_body[row][5] in functional_types:
        functional_types[EAE_PRECLINICAL_body[row][5]] += 1

print(functional_types)


# intronic 399
# 3'5'UTR 180
# downstream 17
# upstream 3
# lncRNA 0
# exonic 144




########################
##### EAE CLINICAL #####
########################

EAE_CLINICAL_total_edited_sites = len(EAE_CLINICAL_body)

print(EAE_CLINICAL_total_edited_sites)
# 546

functional_types = {}

for row in range(0, len(EAE_CLINICAL_body)):
    if EAE_CLINICAL_body[row][5] != '' and EAE_CLINICAL_body[row][5] != 'NA' and EAE_CLINICAL_body[row][5] not in functional_types:
        functional_types[EAE_CLINICAL_body[row][5]] = 1
    elif EAE_CLINICAL_body[row][5] in functional_types:
        functional_types[EAE_CLINICAL_body[row][5]] += 1

print(functional_types)


# intronic 177
# 3'5'UTR 222
# downstream 19
# upstream 4
# lncRNA 0
# exonic 124




######################################################################################
##### WMvsNAIVE SITES FROM COMMON GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
######################################################################################

WMvsNAIVE_total_edited_sites = len(WMvsNAIVE_body)

print(WMvsNAIVE_total_edited_sites)
# 311

functional_types = {}

for row in range(0, len(WMvsNAIVE_body)):
    if WMvsNAIVE_body[row][7] != '' and WMvsNAIVE_body[row][7] != 'NA' and WMvsNAIVE_body[row][7] not in functional_types:
        functional_types[WMvsNAIVE_body[row][7]] = 1
    elif WMvsNAIVE_body[row][7] in functional_types:
        functional_types[WMvsNAIVE_body[row][7]] += 1

print(functional_types)


# intronic 259
# 3'5'UTR 48
# downstream 0
# upstream 3
# lncRNA 1
# exonic 0




##############################################################################################
##### NAWMvsPRECLINICAL SITES FROM COMMON GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
##############################################################################################

NAWMvsPRECLINICAL_total_edited_sites = len(NAWMvsPRECLINICAL_body)

print(NAWMvsPRECLINICAL_total_edited_sites)
# 238

functional_types = {}

for row in range(0, len(NAWMvsPRECLINICAL_body)):
    if NAWMvsPRECLINICAL_body[row][7] != '' and NAWMvsPRECLINICAL_body[row][7] != 'NA' and NAWMvsPRECLINICAL_body[row][7] not in functional_types:
        functional_types[NAWMvsPRECLINICAL_body[row][7]] = 1
    elif NAWMvsPRECLINICAL_body[row][7] in functional_types:
        functional_types[NAWMvsPRECLINICAL_body[row][7]] += 1

print(functional_types)


# intronic 169
# 3'5'UTR 62
# downstream 0
# upstream 0
# lncRNA 7
# exonic 0




#########################################################################################
##### ALvsCLINICAL SITES FROM COMMON GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#########################################################################################

ALvsCLINICAL_total_edited_sites = len(ALvsCLINICAL_body)

print(ALvsCLINICAL_total_edited_sites)
# 87

functional_types = {}

for row in range(0, len(ALvsCLINICAL_body)):
    if ALvsCLINICAL_body[row][7] != '' and ALvsCLINICAL_body[row][7] != 'NA' and ALvsCLINICAL_body[row][7] not in functional_types:
        functional_types[ALvsCLINICAL_body[row][7]] = 1
    elif ALvsCLINICAL_body[row][7] in functional_types:
        functional_types[ALvsCLINICAL_body[row][7]] += 1

print(functional_types)


# intronic 32
# 3'5'UTR 54
# downstream 1
# upstream 0
# lncRNA 0
# exonic 0




#########################################################################################
##### CAvsCLINICAL SITES FROM COMMON GENES AFTER QC MICROGLIA AND SAMPLE FILTRATION #####
#########################################################################################

CAvsCLINICAL_total_edited_sites = len(CAvsCLINICAL_body)

print(CAvsCLINICAL_total_edited_sites)
# 65

functional_types = {}

for row in range(0, len(CAvsCLINICAL_body)):
    if CAvsCLINICAL_body[row][7] != '' and CAvsCLINICAL_body[row][7] != 'NA' and CAvsCLINICAL_body[row][7] not in functional_types:
        functional_types[CAvsCLINICAL_body[row][7]] = 1
    elif CAvsCLINICAL_body[row][7] in functional_types:
        functional_types[CAvsCLINICAL_body[row][7]] += 1

print(functional_types)


# intronic 34
# 3'5'UTR 31
# downstream 0
# upstream 0
# lncRNA 0
# exonic 0




































