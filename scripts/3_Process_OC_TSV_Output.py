#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess

########################################################
### Exclude samples with < 20M uniquely mapped reads ###
###       Curated list after manual review           ###
########################################################

TOTAL_WM_AFTER_QC = 24
TOTAL_NAWM_AFTER_QC = 19
TOTAL_AL_AFTER_QC = 14
TOTAL_CA_AFTER_QC = 11



OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"


##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/3_OC_ANNOTATED"
subprocess.run(cmd, shell=True)


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/2_OC_OUTPUT/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT.vcf.gz.variant.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)


header_line_index = 0

for row_index in range(0, len(vfile)):
    if vfile[row_index][0] == "uid":
        header_line_index = row_index
        break
    else:
        next


header = vfile[header_line_index]

body = vfile[(header_line_index+1):]

# Set the columns to remove
columns_to_remove = ['uid', "note", "all_mappings", "tags", "clinpred.rankscore", "clinvar.disease_refs", "clinvar.rev_stat", "clinvar.id", "encode_tfbs.all", "extra_vcf_info.pos", "extra_vcf_info.ref", "extra_vcf_info.alt", "original_input.chrom", "original_input.pos", "original_input.ref_base", "original_input.alt_base", "hpo.all", "ncbigene.entrez", "vcfinfo.phred", "vcfinfo.filter", "vcfinfo.zygosity", "vcfinfo.alt_reads", "vcfinfo.tot_reads", "vcfinfo.af", "vcfinfo.hap_block", "vcfinfo.hap_strand", 'gnomad3.af_afr', 'gnomad3.af_asj', 'gnomad3.af_eas', 'gnomad3.af_fin', 'gnomad3.af_lat', 'gnomad3.af_nfe', 'gnomad3.af_oth', 'gnomad3.af_sas']

# Retrive indices of columns to remove
columns_to_remove_indices = []

for item in columns_to_remove:
    columns_to_remove_indices.append(header.index(item))

# Remove not needed colums from the body
for row in range(0, len(body)):
    for index in sorted(columns_to_remove_indices, reverse=True):
        del body[row][index]

# Remove not needed columns from the header
for item in columns_to_remove:
    header.remove(item)

# Rename header columns
header[header.index("chrom")] = "Chrom"
header[header.index("pos")] = "hg38_Position"
header[header.index("ref_base")] = "Ref_Base"
header[header.index("alt_base")] = "Alt_Base"
header[header.index("coding")] = "Coding"
header[header.index("hugo")] = "Gene"
header[header.index("transcript")] = "Transcript"
header[header.index("so")] = "Sequence_Ontology"
header[header.index("cchange")] = "cDNA_change"
header[header.index("achange")] = "Protein_Change"
header[header.index("numsample")] = "Samples_Count"
header[header.index("samples")] = "Samples"
header[header.index("hg19.chrom")] = "hg19_Chrom"
header[header.index("hg19.pos")] = "hg19_Pos"
header[header.index("clinpred.score")] = "Clinpred_Score"
header[header.index("clinvar.sig")] = "Clinvar_Significance"
header[header.index("clinvar.disease_names")] = "Clinvar_Disease"
header[header.index("encode_tfbs.factor")] = "Encode_TFBS_Factor"
header[header.index("ensembl_regulatory_build.region")] = "Ensemble_Regulatory_Features"
header[header.index("ensembl_regulatory_build.ensr")] = "Ensemble_Regulatory_Features_ID"
header[header.index("hpo.id")] = "Human_Phenotype_ID"
header[header.index("hpo.term")] = "Human_Phenotype_Term"
header[header.index("ncbigene.ncbi_desc")] = "NCBI_Gene_Description"
header[header.index("gnomad3.af")] = "gnomAD3_AF"
header[header.index("pseudogene.hugo")] = "Pseudogene_Hugo"
header[header.index("pseudogene.transcript")] = "Pseudogene_Transcript"
header[header.index("repeat.repeatclass")] = "Repeat_Class"
header[header.index("repeat.repeatfamily")] = "Repeat_Family"
header[header.index("repeat.repeatname")] = "Repeat_Name"
header[header.index("dbsnp.rsid")] = "dbSNP_rsID"
header[header.index("mirbase.transcript")] = "miRBase_Transcript"
header[header.index("mirbase.id")] = "miRBase_ID"
header[header.index("mirbase.name")] = "miRBase_Name"
header[header.index("mirbase.derives_from")] = "miRBase_DerivesFrom"
header[header.index("ncrna.ncrnaclass")] = "ncRNA_Class"
header[header.index("ncrna.ncrnaname")] = "ncRNA_Name"


# Rename extra_vcf_info columns
for column in range(0, len(header)):
    if header[column].startswith("extra_vcf_info."):
        header[column] = header[column][len("extra_vcf_info."):]

# Get the index of "TYPE" column
type_index = header.index("TYPE")

# Change "TYPE", "STRAND", "MEDIATOR", "DJ", "Homopolymer", "Entropy", "RNAEDIT", "EDTYPE" and "EDREGION" column position in header
# Move "TYPE" column
header[type_index-1], header[type_index] = header[type_index], header[type_index-1]
header[type_index-2], header[type_index-1] = header[type_index-1], header[type_index-2]
header[type_index-3], header[type_index-2] = header[type_index-2], header[type_index-3]
header[type_index-4], header[type_index-3] = header[type_index-3], header[type_index-4]
# Move "STRAND" column
header[type_index], header[type_index+1] = header[type_index+1], header[type_index]
header[type_index-1], header[type_index] = header[type_index], header[type_index-1]
header[type_index-2], header[type_index-1] = header[type_index-1], header[type_index-2]
header[type_index-3], header[type_index-2] = header[type_index-2], header[type_index-3]
# Move "MEDIATOR" column
header[type_index+1], header[type_index+2] = header[type_index+2], header[type_index+1]
header[type_index], header[type_index+1] = header[type_index+1], header[type_index]
header[type_index-1], header[type_index] = header[type_index], header[type_index-1]
header[type_index-2], header[type_index-1] = header[type_index-1], header[type_index-2]
# Move "DJ" column
header[type_index+2], header[type_index+3] = header[type_index+3], header[type_index+2]
header[type_index+1], header[type_index+2] = header[type_index+2], header[type_index+1]
header[type_index], header[type_index+1] = header[type_index+1], header[type_index]
header[type_index-1], header[type_index] = header[type_index], header[type_index-1]
# Move "Homopolymer" column
header[type_index+3], header[type_index+4] = header[type_index+4], header[type_index+3]
header[type_index+2], header[type_index+3] = header[type_index+3], header[type_index+2]
header[type_index+1], header[type_index+2] = header[type_index+2], header[type_index+1]
header[type_index], header[type_index+1] = header[type_index+1], header[type_index]
# Move "Entropy" column
header[type_index+4], header[type_index+5] = header[type_index+5], header[type_index+4]
header[type_index+3], header[type_index+4] = header[type_index+4], header[type_index+3]
header[type_index+2], header[type_index+3] = header[type_index+3], header[type_index+2]
header[type_index+1], header[type_index+2] = header[type_index+2], header[type_index+1]
# Move "RNAEDIT" column
header[type_index+5], header[type_index+6] = header[type_index+6], header[type_index+5]
header[type_index+4], header[type_index+5] = header[type_index+5], header[type_index+4]
header[type_index+3], header[type_index+4] = header[type_index+4], header[type_index+3]
header[type_index+2], header[type_index+3] = header[type_index+3], header[type_index+2]
# Move "EDTYPE" column
header[type_index+6], header[type_index+7] = header[type_index+7], header[type_index+6]
header[type_index+5], header[type_index+6] = header[type_index+6], header[type_index+5]
header[type_index+4], header[type_index+5] = header[type_index+5], header[type_index+4]
header[type_index+3], header[type_index+4] = header[type_index+4], header[type_index+3]
# Move "EDREGION" column
header[type_index+7], header[type_index+8] = header[type_index+8], header[type_index+7]
header[type_index+6], header[type_index+7] = header[type_index+7], header[type_index+6]
header[type_index+5], header[type_index+6] = header[type_index+6], header[type_index+5]
header[type_index+4], header[type_index+5] = header[type_index+5], header[type_index+4]

# Change "TYPE", "STRAND", "MEDIATOR", "DJ", "Homopolymer", "Entropy", "RNAEDIT", "EDTYPE" and "EDREGION" column position in body
for row in range(0, len(body)):
    # Move "TYPE" column
    body[row][type_index-1], body[row][type_index] = body[row][type_index], body[row][type_index-1]
    body[row][type_index-2], body[row][type_index-1] = body[row][type_index-1], body[row][type_index-2]
    body[row][type_index-3], body[row][type_index-2] = body[row][type_index-2], body[row][type_index-3]
    body[row][type_index-4], body[row][type_index-3] = body[row][type_index-3], body[row][type_index-4]
    # Move "STRAND" column
    body[row][type_index], body[row][type_index+1] = body[row][type_index+1], body[row][type_index]
    body[row][type_index-1], body[row][type_index] = body[row][type_index], body[row][type_index-1]
    body[row][type_index-2], body[row][type_index-1] = body[row][type_index-1], body[row][type_index-2]
    body[row][type_index-3], body[row][type_index-2] = body[row][type_index-2], body[row][type_index-3]
    # Move "MEDIATOR" column
    body[row][type_index+1], body[row][type_index+2] = body[row][type_index+2], body[row][type_index+1]
    body[row][type_index], body[row][type_index+1] = body[row][type_index+1], body[row][type_index]
    body[row][type_index-1], body[row][type_index] = body[row][type_index], body[row][type_index-1]
    body[row][type_index-2], body[row][type_index-1] = body[row][type_index-1], body[row][type_index-2]
    # Move "DJ" column
    body[row][type_index+2], body[row][type_index+3] = body[row][type_index+3], body[row][type_index+2]
    body[row][type_index+1], body[row][type_index+2] = body[row][type_index+2], body[row][type_index+1]
    body[row][type_index], body[row][type_index+1] = body[row][type_index+1], body[row][type_index]
    body[row][type_index-1], body[row][type_index] = body[row][type_index], body[row][type_index-1]
    # Move "Homopolymer" column
    body[row][type_index+3], body[row][type_index+4] = body[row][type_index+4], body[row][type_index+3]
    body[row][type_index+2], body[row][type_index+3] = body[row][type_index+3], body[row][type_index+2]
    body[row][type_index+1], body[row][type_index+2] = body[row][type_index+2], body[row][type_index+1]
    body[row][type_index], body[row][type_index+1] = body[row][type_index+1], body[row][type_index]
    # Move "Entropy" column
    body[row][type_index+4], body[row][type_index+5] = body[row][type_index+5], body[row][type_index+4]
    body[row][type_index+3], body[row][type_index+4] = body[row][type_index+4], body[row][type_index+3]
    body[row][type_index+2], body[row][type_index+3] = body[row][type_index+3], body[row][type_index+2]
    body[row][type_index+1], body[row][type_index+2] = body[row][type_index+2], body[row][type_index+1]
    # Move "RNAEDIT" column
    body[row][type_index+5], body[row][type_index+6] = body[row][type_index+6], body[row][type_index+5]
    body[row][type_index+4], body[row][type_index+5] = body[row][type_index+5], body[row][type_index+4]
    body[row][type_index+3], body[row][type_index+4] = body[row][type_index+4], body[row][type_index+3]
    body[row][type_index+2], body[row][type_index+3] = body[row][type_index+3], body[row][type_index+2]
    # Move "EDTYPE" column
    body[row][type_index+6], body[row][type_index+7] = body[row][type_index+7], body[row][type_index+6]
    body[row][type_index+5], body[row][type_index+6] = body[row][type_index+6], body[row][type_index+5]
    body[row][type_index+4], body[row][type_index+5] = body[row][type_index+5], body[row][type_index+4]
    body[row][type_index+3], body[row][type_index+4] = body[row][type_index+4], body[row][type_index+3]
    # Move "EDREGION" column
    body[row][type_index+7], body[row][type_index+8] = body[row][type_index+8], body[row][type_index+7]
    body[row][type_index+6], body[row][type_index+7] = body[row][type_index+7], body[row][type_index+6]
    body[row][type_index+5], body[row][type_index+6] = body[row][type_index+6], body[row][type_index+5]
    body[row][type_index+4], body[row][type_index+5] = body[row][type_index+5], body[row][type_index+4]


# Add extra header columns that show the number or samples in each category
header.insert(10, "WM_Count")
header.insert(11, "NAWM_Count")
header.insert(12, "AL_Count")
header.insert(13, "CA_Count")


# Sample lists
WM_SAMPLES = ["C3_1", "C5_8", "C5_10", "C5_13", "C4_16", "C2_18", "C4_19", "C1_24", "C5_25", "C3_26", "C4_28", "C2_32", "C5_33", "C4_34", "C3_53", "C4_54", "C3_55", "C1_56", "C1_57", "C2_58", "C1_59", "C1_60", "C2_61", "C2_62", "C3_63"]
NAWM_SAMPLES = ["MS3_11", "MS5_15", "MS8_20", "MS1_23", "MS3_27", "MS8_30", "MS3_31", "MS5_35", "MS3_41", "MS4_42", "MS1_43", "MS3_44", "MS1_45", "MS8_46", "MS6_47", "MS9_48", "MS8_49", "MS5_50", "MS6_51", "MS4_52", "MS5_94"]
AL_SAMPLES = ["MS5_64", "MS4_65", "MS6_66", "MS6_67", "MS5_68", "MS1_71", "MS6_72", "MS10_73", "MS1_74", "MS1_76", "MS2_77", "MS10_79", "MS7_86", "MS6_89", "MS4_90", "MS7_91"]
CA_SAMPLES = ["MS7_2", "MS6_3", "MS4_4", "MS8_5", "MS8_6", "MS8_9", "MS3_12", "MS4_14", "MS3_17", "MS4_21", "MS4_22", "MS4_29", "MS3_36", "MS8_37", "MS6_38", "MS10_39", "MS8_40"]


# Add extra body columns with the number or samples in each category
for row in range(0, len(body)):
    NAWM_Count = 0
    WM_Count = 0
    AL_Count = 0
    CA_Count = 0
    row_samples = body[row][11].split(";")
    for sample in row_samples:
        if sample in NAWM_SAMPLES:
            NAWM_Count += 1
        elif sample in WM_SAMPLES:
            WM_Count += 1
        elif sample in AL_SAMPLES:
            AL_Count += 1
        elif sample in CA_SAMPLES:
            CA_Count += 1
    body[row].insert(10, str(eval("WM_Count")))
    body[row].insert(11, str(eval("NAWM_Count")))
    body[row].insert(12, str(eval("AL_Count")))
    body[row].insert(13, str(eval("CA_Count")))


first_sample_column = (header.index("EDREGION")+1)
last_sample_column = (header.index("Human_Phenotype_ID")-1)
# Add NA to empty columns to avoid text wraping in excel
for row in range(0, len(body)):
    for cell in range(0, len(body[row])):
        if (body[row][cell] == '' or body[row][cell] == 'NA') and cell >= first_sample_column and cell <= last_sample_column:
            body[row][cell] = '0'
        # "Coding" column
        elif (body[row][cell] == '' or body[row][cell] == 'NA') and cell == header.index("Coding"):
            body[row][cell] = 'No'
        # Every other empty column
        elif (body[row][cell] == '' or body[row][cell] == 'NA'):
            body[row][cell] = 'NA'


header.insert(14, "WM_Percentage")
header.insert(15, "NAWM_Percentage")
header.insert(16, "AL_Percentage")
header.insert(17, "CA_Percentage")

WM_COUNT_INDEX = header.index("WM_Count")
NAWM_COUNT_INDEX = header.index("NAWM_Count")
AL_COUNT_INDEX = header.index("AL_Count")
CA_COUNT_INDEX = header.index("CA_Count")

for row in range(0, len(body)):
    WM_PERCENTAGE = round((int(body[row][WM_COUNT_INDEX]) / TOTAL_WM_AFTER_QC), 3)
    NAWM_PERCENTAGE = round((int(body[row][NAWM_COUNT_INDEX]) / TOTAL_NAWM_AFTER_QC), 3)
    AL_PERCENTAGE = round((int(body[row][AL_COUNT_INDEX]) / TOTAL_AL_AFTER_QC), 3)
    CA_PERCENTAGE = round((int(body[row][CA_COUNT_INDEX]) / TOTAL_CA_AFTER_QC), 3)
    body[row].insert(14, str(WM_PERCENTAGE))
    body[row].insert(15, str(NAWM_PERCENTAGE))
    body[row].insert(16, str(AL_PERCENTAGE))
    body[row].insert(17, str(CA_PERCENTAGE))



# REFDEPTH Sample lists
AL_REFDEPTH_SAMPLES = ["MS5_64_REFDEPTH", "MS4_65_REFDEPTH", "MS6_66_REFDEPTH", "MS6_67_REFDEPTH", "MS5_68_REFDEPTH", "MS1_71_REFDEPTH", "MS6_72_REFDEPTH", "MS10_73_REFDEPTH", "MS1_74_REFDEPTH", "MS1_76_REFDEPTH", "MS2_77_REFDEPTH", "MS10_79_REFDEPTH", "MS7_86_REFDEPTH", "MS6_89_REFDEPTH", "MS4_90_REFDEPTH", "MS7_91_REFDEPTH"]
NAWM_REFDEPTH_SAMPLES = ["MS3_11_REFDEPTH", "MS5_15_REFDEPTH", "MS8_20_REFDEPTH", "MS1_23_REFDEPTH", "MS3_27_REFDEPTH", "MS8_30_REFDEPTH", "MS3_31_REFDEPTH", "MS5_35_REFDEPTH", "MS3_41_REFDEPTH", "MS4_42_REFDEPTH", "MS1_43_REFDEPTH", "MS3_44_REFDEPTH", "MS1_45_REFDEPTH", "MS8_46_REFDEPTH", "MS6_47_REFDEPTH", "MS9_48_REFDEPTH", "MS8_49_REFDEPTH", "MS5_50_REFDEPTH", "MS6_51_REFDEPTH", "MS4_52_REFDEPTH", "MS5_94_REFDEPTH"]
CA_REFDEPTH_SAMPLES = ["MS7_2_REFDEPTH", "MS6_3_REFDEPTH", "MS4_4_REFDEPTH", "MS8_5_REFDEPTH", "MS8_6_REFDEPTH", "MS8_9_REFDEPTH", "MS3_12_REFDEPTH", "MS4_14_REFDEPTH", "MS3_17_REFDEPTH", "MS4_21_REFDEPTH", "MS4_22_REFDEPTH", "MS4_29_REFDEPTH", "MS3_36_REFDEPTH", "MS8_37_REFDEPTH", "MS6_38_REFDEPTH", "MS10_39_REFDEPTH", "MS8_40_REFDEPTH"]
WM_REFDEPTH_SAMPLES = ["C3_1_REFDEPTH", "C5_8_REFDEPTH", "C5_10_REFDEPTH", "C5_13_REFDEPTH", "C4_16_REFDEPTH", "C2_18_REFDEPTH", "C4_19_REFDEPTH", "C1_24_REFDEPTH", "C5_25_REFDEPTH", "C3_26_REFDEPTH", "C4_28_REFDEPTH", "C2_32_REFDEPTH", "C5_33_REFDEPTH", "C4_34_REFDEPTH", "C3_53_REFDEPTH", "C4_54_REFDEPTH", "C3_55_REFDEPTH", "C1_56_REFDEPTH", "C1_57_REFDEPTH", "C2_58_REFDEPTH", "C1_59_REFDEPTH", "C1_60_REFDEPTH", "C2_61_REFDEPTH", "C2_62_REFDEPTH", "C3_63_REFDEPTH"]

# Initiate sample category indices
WM_REFDEPTH_INDICES = []
CA_REFDEPTH_INDICES = []
NAWM_REFDEPTH_INDICES = []
AL_REFDEPTH_INDICES = []

# Find sample indices 
for cell in range(0, len(header)):
    if header[cell] in WM_REFDEPTH_SAMPLES:
        WM_REFDEPTH_INDICES.append(int(cell))
    elif header[cell] in CA_REFDEPTH_SAMPLES:
        CA_REFDEPTH_INDICES.append(int(cell))
    elif header[cell] in NAWM_REFDEPTH_SAMPLES:
        NAWM_REFDEPTH_INDICES.append(int(cell))
    elif header[cell] in AL_REFDEPTH_SAMPLES:
        AL_REFDEPTH_INDICES.append(int(cell))


# Add extra REFDEPTH values to body file
for row in range(0, len(body)):
    WM_REFDEPTH_SUM = 0
    CA_REFDEPTH_SUM = 0
    NAWM_REFDEPTH_SUM = 0
    AL_REFDEPTH_SUM = 0
    # Change these lines according to samples of interest
    for WM_REFDEPTH_INDEX in WM_REFDEPTH_INDICES:
        WM_REFDEPTH_SUM += int(body[row][WM_REFDEPTH_INDEX])
    for NAWM_REFDEPTH_INDEX in NAWM_REFDEPTH_INDICES:
        NAWM_REFDEPTH_SUM += int(body[row][NAWM_REFDEPTH_INDEX])
    for AL_REFDEPTH_INDEX in AL_REFDEPTH_INDICES:
        AL_REFDEPTH_SUM += int(body[row][AL_REFDEPTH_INDEX])
    for CA_REFDEPTH_INDEX in CA_REFDEPTH_INDICES:
        CA_REFDEPTH_SUM += int(body[row][CA_REFDEPTH_INDEX])
    body[row].insert(18, str(eval("WM_REFDEPTH_SUM")))
    body[row].insert(19, str(eval("NAWM_REFDEPTH_SUM")))
    body[row].insert(20, str(eval("AL_REFDEPTH_SUM")))
    body[row].insert(21, str(eval("CA_REFDEPTH_SUM")))

# Add extra REFDEPTH columns to header file
header.insert(18, "WM_REFDEPTH") 
header.insert(19, "NAWM_REFDEPTH")
header.insert(20, "AL_REFDEPTH")  
header.insert(21, "CA_REFDEPTH")  





 


# ALTDEPTH Sample lists
AL_ALTDEPTH_SAMPLES = ["MS5_64_ALTDEPTH", "MS4_65_ALTDEPTH", "MS6_66_ALTDEPTH", "MS6_67_ALTDEPTH", "MS5_68_ALTDEPTH", "MS1_71_ALTDEPTH", "MS6_72_ALTDEPTH", "MS10_73_ALTDEPTH", "MS1_74_ALTDEPTH", "MS1_76_ALTDEPTH", "MS2_77_ALTDEPTH", "MS10_79_ALTDEPTH", "MS7_86_ALTDEPTH", "MS6_89_ALTDEPTH", "MS4_90_ALTDEPTH", "MS7_91_ALTDEPTH"]
NAWM_ALTDEPTH_SAMPLES = ["MS3_11_ALTDEPTH", "MS5_15_ALTDEPTH", "MS8_20_ALTDEPTH", "MS1_23_ALTDEPTH", "MS3_27_ALTDEPTH", "MS8_30_ALTDEPTH", "MS3_31_ALTDEPTH", "MS5_35_ALTDEPTH", "MS3_41_ALTDEPTH", "MS4_42_ALTDEPTH", "MS1_43_ALTDEPTH", "MS3_44_ALTDEPTH", "MS1_45_ALTDEPTH", "MS8_46_ALTDEPTH", "MS6_47_ALTDEPTH", "MS9_48_ALTDEPTH", "MS8_49_ALTDEPTH", "MS5_50_ALTDEPTH", "MS6_51_ALTDEPTH", "MS4_52_ALTDEPTH", "MS5_94_ALTDEPTH"]
CA_ALTDEPTH_SAMPLES = ["MS7_2_ALTDEPTH", "MS6_3_ALTDEPTH", "MS4_4_ALTDEPTH", "MS8_5_ALTDEPTH", "MS8_6_ALTDEPTH", "MS8_9_ALTDEPTH", "MS3_12_ALTDEPTH", "MS4_14_ALTDEPTH", "MS3_17_ALTDEPTH", "MS4_21_ALTDEPTH", "MS4_22_ALTDEPTH", "MS4_29_ALTDEPTH", "MS3_36_ALTDEPTH", "MS8_37_ALTDEPTH", "MS6_38_ALTDEPTH", "MS10_39_ALTDEPTH", "MS8_40_ALTDEPTH"]
WM_ALTDEPTH_SAMPLES = ["C3_1_ALTDEPTH", "C5_8_ALTDEPTH", "C5_10_ALTDEPTH", "C5_13_ALTDEPTH", "C4_16_ALTDEPTH", "C2_18_ALTDEPTH", "C4_19_ALTDEPTH", "C1_24_ALTDEPTH", "C5_25_ALTDEPTH", "C3_26_ALTDEPTH", "C4_28_ALTDEPTH", "C2_32_ALTDEPTH", "C5_33_ALTDEPTH", "C4_34_ALTDEPTH", "C3_53_ALTDEPTH", "C4_54_ALTDEPTH", "C3_55_ALTDEPTH", "C1_56_ALTDEPTH", "C1_57_ALTDEPTH", "C2_58_ALTDEPTH", "C1_59_ALTDEPTH", "C1_60_ALTDEPTH", "C2_61_ALTDEPTH", "C2_62_ALTDEPTH", "C3_63_ALTDEPTH"]

# Initiate sample category indices
WM_ALTDEPTH_INDICES = []
CA_ALTDEPTH_INDICES = []
NAWM_ALTDEPTH_INDICES = []
AL_ALTDEPTH_INDICES = []

# Find sample indices 
for cell in range(0, len(header)):
    if header[cell] in WM_ALTDEPTH_SAMPLES:
        WM_ALTDEPTH_INDICES.append(int(cell))
    elif header[cell] in CA_ALTDEPTH_SAMPLES:
        CA_ALTDEPTH_INDICES.append(int(cell))
    elif header[cell] in NAWM_ALTDEPTH_SAMPLES:
        NAWM_ALTDEPTH_INDICES.append(int(cell))
    elif header[cell] in AL_ALTDEPTH_SAMPLES:
        AL_ALTDEPTH_INDICES.append(int(cell))


# Add extra ALTDEPTH values to body file
for row in range(0, len(body)):
    WM_ALTDEPTH_SUM = 0
    CA_ALTDEPTH_SUM = 0
    NAWM_ALTDEPTH_SUM = 0
    AL_ALTDEPTH_SUM = 0
    # Change these lines according to samples of interest
    for WM_ALTDEPTH_INDEX in WM_ALTDEPTH_INDICES:
        WM_ALTDEPTH_SUM += int(body[row][WM_ALTDEPTH_INDEX])
    for NAWM_ALTDEPTH_INDEX in NAWM_ALTDEPTH_INDICES:
        NAWM_ALTDEPTH_SUM += int(body[row][NAWM_ALTDEPTH_INDEX])
    for AL_ALTDEPTH_INDEX in AL_ALTDEPTH_INDICES:
        AL_ALTDEPTH_SUM += int(body[row][AL_ALTDEPTH_INDEX])
    for CA_ALTDEPTH_INDEX in CA_ALTDEPTH_INDICES:
        CA_ALTDEPTH_SUM += int(body[row][CA_ALTDEPTH_INDEX])
    body[row].insert(22, str(eval("WM_ALTDEPTH_SUM")))
    body[row].insert(23, str(eval("NAWM_ALTDEPTH_SUM")))
    body[row].insert(24, str(eval("AL_ALTDEPTH_SUM")))
    body[row].insert(25, str(eval("CA_ALTDEPTH_SUM")))

# Add extra ALTDEPTH columns to header file
header.insert(22, "WM_ALTDEPTH") 
header.insert(23, "NAWM_ALTDEPTH")
header.insert(24, "AL_ALTDEPTH")  
header.insert(25, "CA_ALTDEPTH") 




# Create COVERAGE Sample lists
AL_COVERAGE_SAMPLES = ["MS5_64_COVERAGE", "MS4_65_COVERAGE", "MS6_66_COVERAGE", "MS6_67_COVERAGE", "MS5_68_COVERAGE", "MS1_71_COVERAGE", "MS6_72_COVERAGE", "MS10_73_COVERAGE", "MS1_74_COVERAGE", "MS1_76_COVERAGE", "MS2_77_COVERAGE", "MS10_79_COVERAGE", "MS7_86_COVERAGE", "MS6_89_COVERAGE", "MS4_90_COVERAGE", "MS7_91_COVERAGE"]
NAWM_COVERAGE_SAMPLES = ["MS3_11_COVERAGE", "MS5_15_COVERAGE", "MS8_20_COVERAGE", "MS1_23_COVERAGE", "MS3_27_COVERAGE", "MS8_30_COVERAGE", "MS3_31_COVERAGE", "MS5_35_COVERAGE", "MS3_41_COVERAGE", "MS4_42_COVERAGE", "MS1_43_COVERAGE", "MS3_44_COVERAGE", "MS1_45_COVERAGE", "MS8_46_COVERAGE", "MS6_47_COVERAGE", "MS9_48_COVERAGE", "MS8_49_COVERAGE", "MS5_50_COVERAGE", "MS6_51_COVERAGE", "MS4_52_COVERAGE", "MS5_94_COVERAGE"]
CA_COVERAGE_SAMPLES = ["MS7_2_COVERAGE", "MS6_3_COVERAGE", "MS4_4_COVERAGE", "MS8_5_COVERAGE", "MS8_6_COVERAGE", "MS8_9_COVERAGE", "MS3_12_COVERAGE", "MS4_14_COVERAGE", "MS3_17_COVERAGE", "MS4_21_COVERAGE", "MS4_22_COVERAGE", "MS4_29_COVERAGE", "MS3_36_COVERAGE", "MS8_37_COVERAGE", "MS6_38_COVERAGE", "MS10_39_COVERAGE", "MS8_40_COVERAGE"]
WM_COVERAGE_SAMPLES = ["C3_1_COVERAGE", "C5_8_COVERAGE", "C5_10_COVERAGE", "C5_13_COVERAGE", "C4_16_COVERAGE", "C2_18_COVERAGE", "C4_19_COVERAGE", "C1_24_COVERAGE", "C5_25_COVERAGE", "C3_26_COVERAGE", "C4_28_COVERAGE", "C2_32_COVERAGE", "C5_33_COVERAGE", "C4_34_COVERAGE", "C3_53_COVERAGE", "C4_54_COVERAGE", "C3_55_COVERAGE", "C1_56_COVERAGE", "C1_57_COVERAGE", "C2_58_COVERAGE", "C1_59_COVERAGE", "C1_60_COVERAGE", "C2_61_COVERAGE", "C2_62_COVERAGE", "C3_63_COVERAGE"]

# Initiate sample category indices
WM_COVERAGE_INDICES = []
CA_COVERAGE_INDICES = []
NAWM_COVERAGE_INDICES = []
AL_COVERAGE_INDICES = []

# Find sample indices 
for cell in range(0, len(header)):
    if header[cell] in WM_COVERAGE_SAMPLES:
        WM_COVERAGE_INDICES.append(int(cell))
    elif header[cell] in CA_COVERAGE_SAMPLES:
        CA_COVERAGE_INDICES.append(int(cell))
    elif header[cell] in NAWM_COVERAGE_SAMPLES:
        NAWM_COVERAGE_INDICES.append(int(cell))
    elif header[cell] in AL_COVERAGE_SAMPLES:
        AL_COVERAGE_INDICES.append(int(cell))


# Add extra COVERAGE values to body file
for row in range(0, len(body)):
    WM_COVERAGE_SUM = 0
    CA_COVERAGE_SUM = 0
    NAWM_COVERAGE_SUM = 0
    AL_COVERAGE_SUM = 0
    # Change these lines according to samples of interest
    for WM_COVERAGE_INDEX in WM_COVERAGE_INDICES:
        WM_COVERAGE_SUM += int(body[row][WM_COVERAGE_INDEX])
    for NAWM_COVERAGE_INDEX in NAWM_COVERAGE_INDICES:
        NAWM_COVERAGE_SUM += int(body[row][NAWM_COVERAGE_INDEX])
    for AL_COVERAGE_INDEX in AL_COVERAGE_INDICES:
        AL_COVERAGE_SUM += int(body[row][AL_COVERAGE_INDEX])
    for CA_COVERAGE_INDEX in CA_COVERAGE_INDICES:
        CA_COVERAGE_SUM += int(body[row][CA_COVERAGE_INDEX])
    body[row].insert(26, str(eval("WM_COVERAGE_SUM")))
    body[row].insert(27, str(eval("NAWM_COVERAGE_SUM")))
    body[row].insert(28, str(eval("AL_COVERAGE_SUM")))
    body[row].insert(29, str(eval("CA_COVERAGE_SUM")))

# Add extra COVERAGE columns to header file
header.insert(26, "WM_COVERAGE") 
header.insert(27, "NAWM_COVERAGE")
header.insert(28, "AL_COVERAGE")  
header.insert(29, "CA_COVERAGE")  




# Create FREQUENCY Sample lists
WM_FREQUENCY_SAMPLES = ["C3_1_FREQUENCY", "C5_8_FREQUENCY", "C5_10_FREQUENCY", "C5_13_FREQUENCY", "C4_16_FREQUENCY", "C2_18_FREQUENCY", "C4_19_FREQUENCY", "C1_24_FREQUENCY", "C5_25_FREQUENCY", "C3_26_FREQUENCY", "C4_28_FREQUENCY", "C2_32_FREQUENCY", "C5_33_FREQUENCY", "C4_34_FREQUENCY", "C3_53_FREQUENCY", "C4_54_FREQUENCY", "C3_55_FREQUENCY", "C1_56_FREQUENCY", "C1_57_FREQUENCY", "C2_58_FREQUENCY", "C1_59_FREQUENCY", "C1_60_FREQUENCY", "C2_61_FREQUENCY", "C2_62_FREQUENCY", "C3_63_FREQUENCY"]
NAWM_FREQUENCY_SAMPLES = ["MS3_11_FREQUENCY", "MS5_15_FREQUENCY", "MS8_20_FREQUENCY", "MS1_23_FREQUENCY", "MS3_27_FREQUENCY", "MS8_30_FREQUENCY", "MS3_31_FREQUENCY", "MS5_35_FREQUENCY", "MS3_41_FREQUENCY", "MS4_42_FREQUENCY", "MS1_43_FREQUENCY", "MS3_44_FREQUENCY", "MS1_45_FREQUENCY", "MS8_46_FREQUENCY", "MS6_47_FREQUENCY", "MS9_48_FREQUENCY", "MS8_49_FREQUENCY", "MS5_50_FREQUENCY", "MS6_51_FREQUENCY", "MS4_52_FREQUENCY", "MS5_94_FREQUENCY"]
AL_FREQUENCY_SAMPLES = ["MS5_64_FREQUENCY", "MS4_65_FREQUENCY", "MS6_66_FREQUENCY", "MS6_67_FREQUENCY", "MS5_68_FREQUENCY", "MS1_71_FREQUENCY", "MS6_72_FREQUENCY", "MS10_73_FREQUENCY", "MS1_74_FREQUENCY", "MS1_76_FREQUENCY", "MS2_77_FREQUENCY", "MS10_79_FREQUENCY", "MS7_86_FREQUENCY", "MS6_89_FREQUENCY", "MS4_90_FREQUENCY", "MS7_91_FREQUENCY"]
CA_FREQUENCY_SAMPLES = ["MS7_2_FREQUENCY", "MS6_3_FREQUENCY", "MS4_4_FREQUENCY", "MS8_5_FREQUENCY", "MS8_6_FREQUENCY", "MS8_9_FREQUENCY", "MS3_12_FREQUENCY", "MS4_14_FREQUENCY", "MS3_17_FREQUENCY", "MS4_21_FREQUENCY", "MS4_22_FREQUENCY", "MS4_29_FREQUENCY", "MS3_36_FREQUENCY", "MS8_37_FREQUENCY", "MS6_38_FREQUENCY", "MS10_39_FREQUENCY", "MS8_40_FREQUENCY"]

# Initiate sample category indices
WM_FREQUENCY_INDICES = []
NAWM_FREQUENCY_INDICES = []
AL_FREQUENCY_INDICES = []
CA_FREQUENCY_INDICES = []

# Find sample indices 
for cell in range(0, len(header)):
    if header[cell] in WM_FREQUENCY_SAMPLES:
        WM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in NAWM_FREQUENCY_SAMPLES:
        NAWM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in AL_FREQUENCY_SAMPLES:
        AL_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in CA_FREQUENCY_SAMPLES:
        CA_FREQUENCY_INDICES.append(int(cell))
    

# Add extra mean frequency values to body file
for row in range(0, len(body)):
    WM_FREQUENCY_SUM = 0
    NAWM_FREQUENCY_SUM = 0
    AL_FREQUENCY_SUM = 0
    CA_FREQUENCY_SUM = 0
    # Change these lines according to samples of interest
    for WM_FREQUENCY_INDEX in WM_FREQUENCY_INDICES:
        WM_FREQUENCY_SUM += float(body[row][WM_FREQUENCY_INDEX])
    for NAWM_FREQUENCY_INDEX in NAWM_FREQUENCY_INDICES:
        NAWM_FREQUENCY_SUM += float(body[row][NAWM_FREQUENCY_INDEX])
    for AL_FREQUENCY_INDEX in AL_FREQUENCY_INDICES:
        AL_FREQUENCY_SUM += float(body[row][AL_FREQUENCY_INDEX])
    for CA_FREQUENCY_INDEX in CA_FREQUENCY_INDICES:
        CA_FREQUENCY_SUM += float(body[row][CA_FREQUENCY_INDEX])
    body[row].insert(30, str(round(eval("WM_FREQUENCY_SUM")/eval("TOTAL_WM_AFTER_QC"), 3)) )
    body[row].insert(31, str(round(eval("NAWM_FREQUENCY_SUM")/eval("TOTAL_NAWM_AFTER_QC"), 3)) )
    body[row].insert(32, str(round(eval("AL_FREQUENCY_SUM")/eval("TOTAL_AL_AFTER_QC"), 3)) )
    body[row].insert(33, str(round(eval("CA_FREQUENCY_SUM")/eval("TOTAL_CA_AFTER_QC"), 3)) )

# Add extra mean frequency columns to header file
header.insert(30, "WM_MEAN_FREQUENCY") 
header.insert(31, "NAWM_MEAN_FREQUENCY")
header.insert(32, "AL_MEAN_FREQUENCY")
header.insert(33, "CA_MEAN_FREQUENCY")  




WM_MEAN_FREQUENCY_INDEX = header.index("WM_MEAN_FREQUENCY")
NAWM_MEAN_FREQUENCY_INDEX = header.index("NAWM_MEAN_FREQUENCY")
AL_MEAN_FREQUENCY_INDEX = header.index("AL_MEAN_FREQUENCY")
CA_MEAN_FREQUENCY_INDEX = header.index("CA_MEAN_FREQUENCY")


# Add extra EDITING PERCENTAGE values to body file
for row in range(0, len(body)):
    if (float(body[row][WM_MEAN_FREQUENCY_INDEX]) > 0) or (float(body[row][AL_MEAN_FREQUENCY_INDEX]) > 0):
        body[row].insert(34, str(round(float(body[row][AL_MEAN_FREQUENCY_INDEX])/(float(body[row][AL_MEAN_FREQUENCY_INDEX])+float(body[row][WM_MEAN_FREQUENCY_INDEX])), 3)) )
    else:
        body[row].insert(34, "NA")
    if (float(body[row][WM_MEAN_FREQUENCY_INDEX]) > 0) or (float(body[row][CA_MEAN_FREQUENCY_INDEX]) > 0):
        body[row].insert(35, str(round(float(body[row][CA_MEAN_FREQUENCY_INDEX])/(float(body[row][CA_MEAN_FREQUENCY_INDEX])+float(body[row][WM_MEAN_FREQUENCY_INDEX])), 3)) )
    else:
        body[row].insert(35, "NA")
    if (float(body[row][NAWM_MEAN_FREQUENCY_INDEX]) > 0) or (float(body[row][AL_MEAN_FREQUENCY_INDEX]) > 0):
        body[row].insert(36, str(round(float(body[row][AL_MEAN_FREQUENCY_INDEX])/(float(body[row][AL_MEAN_FREQUENCY_INDEX])+float(body[row][NAWM_MEAN_FREQUENCY_INDEX])), 3)) )
    else:
        body[row].insert(36, "NA")
    if (float(body[row][NAWM_MEAN_FREQUENCY_INDEX]) > 0) or (float(body[row][CA_MEAN_FREQUENCY_INDEX]) > 0):
        body[row].insert(37, str(round(float(body[row][CA_MEAN_FREQUENCY_INDEX])/(float(body[row][CA_MEAN_FREQUENCY_INDEX])+float(body[row][NAWM_MEAN_FREQUENCY_INDEX])), 3)) )
    else:
        body[row].insert(37, "NA")

# Add extra EDITING PERCENTAGE columns to header file
header.insert(34, "WMvsAL_AL_EDITING_PERCENTAGE") 
header.insert(35, "WMvsCA_CA_EDITING_PERCENTAGE")
header.insert(36, "NAWMvsAL_AL_EDITING_PERCENTAGE")
header.insert(37, "NAWMvsCA_CA_EDITING_PERCENTAGE")





fvfile = []
fvfile.append(header)
for row in range(0, len(body)):
    fvfile.append(body[row])

variant_file.close()

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/3_OC_ANNOTATED/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_Annotated.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()










