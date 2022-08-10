##############################
##############################
### MS GSM4114296 ANALYSIS ###
##############################
##############################

import csv
import subprocess

STUDY_CATEGORY = "MS"
STUDY_ACCESSION = "GSM4114296"

SPRINT_OUTPUT_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/SPRINT/output"
SCRIPTS_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/scripts"
RESOURCES_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/resources"
REFERENCE_GENOME_FILE = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/resources/hg38-noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta"

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/postprocessing"


##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/1_VCF"
subprocess.run(cmd, shell=True)
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/2_OC_OUTPUT"
subprocess.run(cmd, shell=True)

########################################################
### Exclude samples with < 20M uniquely mapped reads ###
###       Curated list after manual review           ###
########################################################
SAMPLES_TO_EXCLUDE_DUE_TO_QC = {"MS4_21":"CA", "MS10_79":"AL", "MS6_38":"CA", "MS3_27":"NAWM", "MS4_29":"CA", "MS4_52":"NAWM", "MS8_37":"CA", "MS2_77":"AL", "MS4_14":"CA", "C4_16":"WM", "MS3_17":"CA"}

TOTAL_WM_INITIAL = 25
TOTAL_NAWM_INITIAL = 21
TOTAL_AL_INITIAL = 16
TOTAL_CA_INITIAL = 17

TOTAL_WM_AFTER_QC = 24
TOTAL_NAWM_AFTER_QC = 19
TOTAL_AL_AFTER_QC = 14
TOTAL_CA_AFTER_QC = 11

###############################################
### SAMPLE NAMES TO INCLUDE IN OUR ANALYSIS ###
###############################################
WM_SAMPLES = ["C3_1", "C5_8", "C5_10", "C5_13", "C4_16", "C2_18", "C4_19", "C1_24", "C5_25", "C3_26", "C4_28", "C2_32", "C5_33", "C4_34", "C3_53", "C4_54", "C3_55", "C1_56", "C1_57", "C2_58", "C1_59", "C1_60", "C2_61", "C2_62", "C3_63"]
NAWM_SAMPLES = ["MS3_11", "MS5_15", "MS8_20", "MS1_23", "MS3_27", "MS8_30", "MS3_31", "MS5_35", "MS3_41", "MS4_42", "MS1_43", "MS3_44", "MS1_45", "MS8_46", "MS6_47", "MS9_48", "MS8_49", "MS5_50", "MS6_51", "MS4_52", "MS5_94"]
AL_SAMPLES = ["MS5_64", "MS4_65", "MS6_66", "MS6_67", "MS5_68", "MS1_71", "MS6_72", "MS10_73", "MS1_74", "MS1_76", "MS2_77", "MS10_79", "MS7_86", "MS6_89", "MS4_90", "MS7_91"]
CA_SAMPLES = ["MS7_2", "MS6_3", "MS4_4", "MS8_5", "MS8_6", "MS8_9", "MS3_12", "MS4_14", 'MS3_17', "MS4_21", 'MS4_22', "MS4_29", "MS3_36", "MS8_37", "MS6_38", "MS10_39", "MS8_40"]

SAMPLE_NAMES = []
for WM_sample in WM_SAMPLES:
    if WM_sample not in SAMPLES_TO_EXCLUDE_DUE_TO_QC:
        SAMPLE_NAMES.append(WM_sample)
for NAWM_sample in NAWM_SAMPLES:
    if NAWM_sample not in SAMPLES_TO_EXCLUDE_DUE_TO_QC:
        SAMPLE_NAMES.append(NAWM_sample)
for AL_sample in AL_SAMPLES:
    if AL_sample not in SAMPLES_TO_EXCLUDE_DUE_TO_QC:
        SAMPLE_NAMES.append(AL_sample)
for CA_sample in CA_SAMPLES:
    if CA_sample not in SAMPLES_TO_EXCLUDE_DUE_TO_QC:
        SAMPLE_NAMES.append(CA_sample)


### KEEP THE LIST ORDERED, FROM NORMAL -> DISEASE TYPE
SAMPLE_TYPES = ["WM", "NAWM", "AL", "CA"]


for SAMPLE_NAME in SAMPLE_NAMES:
    if SAMPLE_NAME in WM_SAMPLES:
        SAMPLE_TYPE = "WM"
    elif SAMPLE_NAME in NAWM_SAMPLES:
        SAMPLE_TYPE = "NAWM"
    elif SAMPLE_NAME in AL_SAMPLES:
        SAMPLE_TYPE = "AL"
    elif SAMPLE_NAME in CA_SAMPLES:
        SAMPLE_TYPE = "CA"

    HYPER_IDENTIFIED_FILE = f"{SPRINT_OUTPUT_PATH}/{SAMPLE_NAME}/SPRINT_identified_hyper.res"
    REGULAR_IDENTIFIED_FILE = f"{SPRINT_OUTPUT_PATH}/{SAMPLE_NAME}/SPRINT_identified_regular.res"

    cmd = f"mkdir -p {SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}"
    subprocess.run(cmd, shell=True)

    fo=open(f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.ALL.SPRINT.tsv",'w')

    output=[]

    r=set()

    fi=open(REGULAR_IDENTIFIED_FILE)

    fi.readline()
    for line in fi:
        seq=line.rstrip().split('\t')
        output.append([seq[0],int(seq[2]),line.rstrip(),'regular'])
        r.add(seq[0]+'|'+seq[2]+'|'+seq[3]) 

    fi.close()


    h=set()

    fi=open(HYPER_IDENTIFIED_FILE)

    fi.readline()
    for line in fi:
        seq=line.rstrip().split('\t')
        output.append([seq[0],int(seq[2]),line.rstrip(),'hyper'])    
        h.add(seq[0]+'|'+seq[2]+'|'+seq[3]) 
    fi.close()

    over = r & h


    output.sort()
    old=set()
    for one in output:
        seq=one[2].split('\t')
        if seq[0]+'|'+seq[2]+'|'+seq[3] not in old:
            old.add(seq[0]+'|'+seq[2]+'|'+seq[3])

            if seq[0]+'|'+seq[2]+'|'+seq[3] in over:
                out='\t'.join(seq[0:4])+'\t'+'regular_and_hyper\t'+seq[5]+'\t'+seq[6]+'\n'
                fo.write(out)
            
            else:
                out='\t'.join(seq[0:4])+'\t'+one[3]+'\t'+seq[5]+'\t'+seq[6]+'\n'
                fo.write(out)




    # Read ALL SPRINT output tables ##
    with open(f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.ALL.SPRINT.tsv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        duplets = {}
        for row in csv_reader:
            if row[3] not in duplets:
                duplets[row[3]] = 1
            else:
                duplets[row[3]] += 1
        

    with open(f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.ALL.SPRINT.duplets.tsv", mode='w') as duplets_file:
        duplets_file = csv.writer(duplets_file, delimiter='\t', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        for duplet, value in duplets.items():
            duplets_file.writerow([duplet, value])



    # Instanciate VEP format dict ##
    vcfdict = {}
    vcfdict['chromosome'] = []
    vcfdict['start'] = []
    vcfdict['end'] = []
    vcfdict['allele'] = []
    vcfdict['coverage'] = []
    vcfdict['basecount'] = []
    vcfdict['frequency'] = []
    vcfdict['type'] = []
    vcfdict['strand'] = []
    vcfdict['mediator'] = []


    hg38_chr_names_to_include = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY"]


    # Read ALL SPRINT output tables ##
    with open(f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.ALL.SPRINT.tsv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0
        for row in csv_reader:
            line_count += 1
            if row[0] == 'Region':
                next
            elif row[0] not in hg38_chr_names_to_include:
                vcfdict['chromosome'].append("NA")
                vcfdict['start'].append("NA")
                vcfdict['end'].append("NA")
                vcfdict['allele'].append("NA")
                vcfdict['basecount'].append("NA")
                vcfdict['coverage'].append("NA")
                vcfdict['frequency'].append("NA")
                vcfdict['type'].append("NA")
                vcfdict['strand'].append("NA")
                vcfdict['mediator'].append("NA")
            else:
                vcfdict['chromosome'].append(row[0])
                vcfdict['start'].append(row[2])
                vcfdict['end'].append(row[2])
                vcfdict['allele'].append(row[3])
                depthrow = row[6].split(":")
                vcfdict['basecount'].append(depthrow[0])
                vcfdict['coverage'].append(depthrow[1])
                vcfdict['frequency'].append(str(round(int(depthrow[0]) / int(depthrow[1]), 3)))
                vcfdict['type'].append(row[4])
                if row[5] == ".":
                    vcfdict['strand'].append("NA")
                else:
                    vcfdict['strand'].append(row[5])
                if row[3] == 'AG' or row[3] == 'TC':
                    vcfdict['mediator'].append("ADAR")
                elif row[3] == 'CT' or row[3] == 'GA':
                    vcfdict['mediator'].append("APOBEC")
                else:
                    vcfdict['mediator'].append("NA")


    # Write COVERAGE AND FREQUENCY FILTERED VCF file ##
    with open(f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.vcf", mode='w') as vcfcandidates_file:
        vcfcandidates_file = csv.writer(vcfcandidates_file, delimiter='\t', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        vcfcandidates_file.writerow(['##fileformat=VCFv4.1'])
        vcfcandidates_file.writerow(['##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set'])
        vcfcandidates_file.writerow(['##INFO=<ID=' + SAMPLE_NAME + '_REFDEPTH,Number=1,Type=Integer,Description="Reference Base Depth">'])
        vcfcandidates_file.writerow(['##INFO=<ID=' + SAMPLE_NAME + '_ALTDEPTH,Number=1,Type=Integer,Description="Alternative Base Depth">'])
        vcfcandidates_file.writerow(['##INFO=<ID=' + SAMPLE_NAME + '_FREQUENCY,Number=1,Type=Float,Description="Allele Frequency">'])
        vcfcandidates_file.writerow(['##INFO=<ID=' + SAMPLE_NAME + '_COVERAGE,Number=1,Type=Integer,Description="Site Coverage">'])
        vcfcandidates_file.writerow(['##INFO=<ID=TYPE,Number=1,Type=String,Description="Regular and/or Hyper RESs">'])
        vcfcandidates_file.writerow(['##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand Oriendation">'])
        vcfcandidates_file.writerow(['##INFO=<ID=MEDIATOR,Number=1,Type=String,Description="Possible RNA Editing Mediator">'])
        vcfcandidates_file.writerow(['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'])
        vcfcandidates_file.writerow(['##contig=<ID=chr1>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr2>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr3>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr4>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr5>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr6>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr7>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr8>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr9>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr10>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr11>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr12>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr13>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr14>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr15>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr16>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr17>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr18>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr19>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr20>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr21>'])
        vcfcandidates_file.writerow(['##contig=<ID=chr22>'])
        vcfcandidates_file.writerow(['##contig=<ID=chrX>'])
        vcfcandidates_file.writerow(['##contig=<ID=chrY>'])
        vcfcandidates_file.writerow(['##contig=<ID=chrM>'])
        vcfcandidates_file.writerow(['##SAMPLE=<ID=' + SAMPLE_NAME + ',TYPE=' + SAMPLE_TYPE + '>'])
        vcfcandidates_file.writerow(['#CHROM', 'POS', 'ID',  'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', SAMPLE_NAME])
        for i in range(0, line_count):
            if vcfdict['chromosome'][i] == "NA":
                next
            ### Require AF >= 10% and Coverage >= 20. This implies ALT Depth >= 2. ###
            elif float(vcfdict['frequency'][i]) >= 0.10 and int(vcfdict['coverage'][i]) >= 20:
                if (vcfdict['allele'][i][0] == "T" and vcfdict['allele'][i][1] == "C") or (vcfdict['allele'][i][0] == "A" and vcfdict['allele'][i][1] == "G") or (vcfdict['allele'][i][0] == "C" and vcfdict['allele'][i][1] == "T") or (vcfdict['allele'][i][0] == "G" and vcfdict['allele'][i][1] == "A"):
                    infofield = 'TYPE=' + vcfdict['type'][i]+ ';STRAND=' + vcfdict['strand'][i] + ';MEDIATOR=' + vcfdict['mediator'][i] + ';' + SAMPLE_NAME + '_REFDEPTH=' + str(int(vcfdict['coverage'][i]) - int(vcfdict['basecount'][i])) + ';' + SAMPLE_NAME + '_ALTDEPTH=' + vcfdict['basecount'][i] + ';' + SAMPLE_NAME + '_FREQUENCY=' + vcfdict['frequency'][i] + ';' + SAMPLE_NAME + '_COVERAGE=' + vcfdict['coverage'][i]
                    vcfcandidates_file.writerow([vcfdict['chromosome'][i], vcfdict['start'][i], '.', vcfdict['allele'][i][0], vcfdict['allele'][i][1], '.', '.', infofield, 'GT', '0/1'])

    

    annotate_DJ_path = f'{SCRIPTS_PATH}/annotate_DJ.py'
    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.vcf"
    splice_bed = f'{RESOURCES_PATH}/SpliceSites/ref_annot.splice_adj.bed.gz'
    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.vcf"
    temp_dir = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/temp"
    cmd = f"python {annotate_DJ_path} --input_vcf {input_vcf} --splice_bed {splice_bed} --output_vcf {output_vcf} --temp_dir {temp_dir}"
    subprocess.run(cmd, shell=True)
    cmd = f"rm -rf {temp_dir}"
    subprocess.run(cmd, shell=True)
    cmd = f"rm -rf {input_vcf}"
    subprocess.run(cmd, shell=True)
    
    annotate_HOMOPOLYMER_path = f'{SCRIPTS_PATH}/annotate_entropy_n_homopolymers.py'
    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.vcf"
    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.vcf"
    cmd = f"python {annotate_HOMOPOLYMER_path} --input_vcf {input_vcf} --ref_genome_fa {REFERENCE_GENOME_FILE} --output_vcf {output_vcf}"
    subprocess.run(cmd, shell=True)
    cmd = f"rm -rf {input_vcf}"
    subprocess.run(cmd, shell=True)

    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.vcf"
    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.sorted.vcf.gz"
    cmd = f"bcftools sort -Oz -o {output_vcf} {input_vcf}"
    subprocess.run(cmd, shell=True)
    cmd = f"rm -rf {input_vcf}"
    subprocess.run(cmd, shell=True)

    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.sorted.vcf.gz"
    cmd = f"bcftools index -f --threads 15 {output_vcf}"
    subprocess.run(cmd, shell=True)

    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.sorted.vcf.gz"
    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.RNADB.sorted.vcf.gz"
    rediportal_db = f"{RESOURCES_PATH}/REDIPortal/REDIPortal_hg38_Aug_2022_Processed.vcf.gz"
    cmd = f"bcftools annotate --output-type z --annotations {rediportal_db} --columns 'INFO/RNAEDIT,INFO/EDTYPE,INFO/EDREGION' --output {output_vcf} {input_vcf}"
    subprocess.run(cmd, shell=True)

    output_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.RNADB.sorted.vcf.gz"
    cmd = f"bcftools index -f --threads 15 {output_vcf}"
    subprocess.run(cmd, shell=True)

    cmd = f"rm -rf {SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.sorted.vcf.gz"
    subprocess.run(cmd, shell=True)
    
    cmd = f"rm -rf {SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.sorted.vcf.gz.csi"
    subprocess.run(cmd, shell=True)

    destination = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{SAMPLE_TYPE}"
    cmd = f"mkdir -p {destination}"
    subprocess.run(cmd, shell=True)

    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.RNADB.sorted.vcf.gz"
    cmd = f"cp {input_vcf} {destination}"
    subprocess.run(cmd, shell=True)

    input_vcf = f"{SPRINT_OUTPUT_PATH}_processed/{SAMPLE_NAME}/{SAMPLE_NAME}.A2IC2U.SPRINT.DJ.HP.RNADB.sorted.vcf.gz.csi"
    cmd = f"cp {input_vcf} {destination}"
    subprocess.run(cmd, shell=True)


for sample_type in SAMPLE_TYPES:
    output_vcf = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{sample_type}/{STUDY_CATEGORY}_{STUDY_ACCESSION}_{sample_type}_A2IC2U_SPRINT.vcf.gz"
    vcf_files = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{sample_type}/*.vcf.gz"
    cmd = f"bcftools merge --missing-to-ref -Oz -o {output_vcf} {vcf_files}"
    subprocess.run(cmd, shell=True)

for sample_type in SAMPLE_TYPES:
    sample_type_file = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{sample_type}/{STUDY_CATEGORY}_{STUDY_ACCESSION}_{sample_type}_A2IC2U_SPRINT.vcf.gz"
    cmd = f"bcftools index -f --threads 15 {sample_type_file}"
    subprocess.run(cmd, shell=True)



################################
### CREATE COMBINED VCF FILE ###
################################

COMBINED_FILE_CONTAINING_ALL_SAMPLES = ""

for sample_type in SAMPLE_TYPES:
    if COMBINED_FILE_CONTAINING_ALL_SAMPLES == "":
        COMBINED_FILE_CONTAINING_ALL_SAMPLES += sample_type
    else:
        COMBINED_FILE_CONTAINING_ALL_SAMPLES += "vs"
        COMBINED_FILE_CONTAINING_ALL_SAMPLES += sample_type


combined_file_directory = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{COMBINED_FILE_CONTAINING_ALL_SAMPLES}/"
cmd = f"mkdir -p {combined_file_directory}"
subprocess.run(cmd, shell=True)

combined_file = f"{OUTPUT_POSTPROCESSING_PATH}/1_VCF/{COMBINED_FILE_CONTAINING_ALL_SAMPLES}/{STUDY_CATEGORY}_{STUDY_ACCESSION}_{COMBINED_FILE_CONTAINING_ALL_SAMPLES}_A2IC2U_SPRINT.vcf.gz"

cmd = f"bcftools merge --missing-to-ref -Oz -o {combined_file}"

for sample_type in SAMPLE_TYPES:
    cmd += f" {OUTPUT_POSTPROCESSING_PATH}/1_VCF/{sample_type}/{STUDY_CATEGORY}_{STUDY_ACCESSION}_{sample_type}_A2IC2U_SPRINT.vcf.gz"

subprocess.run(cmd, shell=True)

cmd = f"bcftools index -f --threads 15 {combined_file}"
subprocess.run(cmd, shell=True)










