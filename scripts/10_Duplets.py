##############################
##############################
### MS GSE138614 ANALYSIS ###
##############################
##############################

import csv
import subprocess
import matplotlib.pyplot as plt


SPRINT_PROCESSED_OUTPUT_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/SPRINT/output_processed"

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"



##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/10_DUPLETS"
subprocess.run(cmd, shell=True)


WM_SAMPLE_NAMES = ["C3_1", "C5_8", "C5_10", "C5_13", "C2_18", "C4_19", "C1_24", "C5_25", "C3_26", "C4_28", "C2_32", "C5_33", "C4_34", "C3_53", "C4_54", "C3_55", "C1_56", "C1_57", "C2_58", "C1_59", "C1_60", "C2_61", "C2_62", "C3_63"]

NAWM_SAMPLE_NAMES = ["MS3_11", "MS5_15", "MS8_20", "MS1_23", "MS8_30", "MS3_31", "MS5_35", "MS3_41", "MS4_42", "MS1_43", "MS3_44", "MS1_45", "MS8_46", "MS6_47", "MS9_48", "MS8_49", "MS5_50", "MS6_51", "MS5_94"]

AL_SAMPLE_NAMES = ["MS5_64", "MS4_65", "MS6_66", "MS6_67", "MS5_68", "MS1_71", "MS6_72", "MS10_73", "MS1_74", "MS1_76", "MS7_86", "MS6_89", "MS4_90", "MS7_91"]

CA_SAMPLE_NAMES = ["MS7_2", "MS6_3", "MS4_4", "MS8_5", "MS8_6", "MS8_9", "MS3_12", 'MS4_22', "MS3_36", "MS10_39", "MS8_40"]

WM_TOTAL = 0
WM_DUPLETS = {}
WM_DUPLETS["TC"] = 0
WM_DUPLETS["AG"] = 0
WM_DUPLETS["GA"] = 0
WM_DUPLETS["CT"] = 0
WM_DUPLETS["AT"] = 0
WM_DUPLETS["GT"] = 0
WM_DUPLETS["TA"] = 0
WM_DUPLETS["CG"] = 0
WM_DUPLETS["CA"] = 0
WM_DUPLETS["AC"] = 0
WM_DUPLETS["TG"] = 0
WM_DUPLETS["GC"] = 0

WM_DUPLETS_PERCENT = {}
WM_DUPLETS_PERCENT["TC"] = 0
WM_DUPLETS_PERCENT["AG"] = 0
WM_DUPLETS_PERCENT["GA"] = 0
WM_DUPLETS_PERCENT["CT"] = 0
WM_DUPLETS_PERCENT["AT"] = 0
WM_DUPLETS_PERCENT["GT"] = 0
WM_DUPLETS_PERCENT["TA"] = 0
WM_DUPLETS_PERCENT["CG"] = 0
WM_DUPLETS_PERCENT["CA"] = 0
WM_DUPLETS_PERCENT["AC"] = 0
WM_DUPLETS_PERCENT["TG"] = 0
WM_DUPLETS_PERCENT["GC"] = 0


NAWM_TOTAL = 0
NAWM_DUPLETS = {}
NAWM_DUPLETS["TC"] = 0
NAWM_DUPLETS["AG"] = 0
NAWM_DUPLETS["GA"] = 0
NAWM_DUPLETS["CT"] = 0
NAWM_DUPLETS["AT"] = 0
NAWM_DUPLETS["GT"] = 0
NAWM_DUPLETS["TA"] = 0
NAWM_DUPLETS["CG"] = 0
NAWM_DUPLETS["CA"] = 0
NAWM_DUPLETS["AC"] = 0
NAWM_DUPLETS["TG"] = 0
NAWM_DUPLETS["GC"] = 0

NAWM_DUPLETS_PERCENT = {}
NAWM_DUPLETS_PERCENT["TC"] = 0
NAWM_DUPLETS_PERCENT["AG"] = 0
NAWM_DUPLETS_PERCENT["GA"] = 0
NAWM_DUPLETS_PERCENT["CT"] = 0
NAWM_DUPLETS_PERCENT["AT"] = 0
NAWM_DUPLETS_PERCENT["GT"] = 0
NAWM_DUPLETS_PERCENT["TA"] = 0
NAWM_DUPLETS_PERCENT["CG"] = 0
NAWM_DUPLETS_PERCENT["CA"] = 0
NAWM_DUPLETS_PERCENT["AC"] = 0
NAWM_DUPLETS_PERCENT["TG"] = 0
NAWM_DUPLETS_PERCENT["GC"] = 0


AL_TOTAL = 0
AL_DUPLETS = {}
AL_DUPLETS["TC"] = 0
AL_DUPLETS["AG"] = 0
AL_DUPLETS["GA"] = 0
AL_DUPLETS["CT"] = 0
AL_DUPLETS["AT"] = 0
AL_DUPLETS["GT"] = 0
AL_DUPLETS["TA"] = 0
AL_DUPLETS["CG"] = 0
AL_DUPLETS["CA"] = 0
AL_DUPLETS["AC"] = 0
AL_DUPLETS["TG"] = 0
AL_DUPLETS["GC"] = 0

AL_DUPLETS_PERCENT = {}
AL_DUPLETS_PERCENT["TC"] = 0
AL_DUPLETS_PERCENT["AG"] = 0
AL_DUPLETS_PERCENT["GA"] = 0
AL_DUPLETS_PERCENT["CT"] = 0
AL_DUPLETS_PERCENT["AT"] = 0
AL_DUPLETS_PERCENT["GT"] = 0
AL_DUPLETS_PERCENT["TA"] = 0
AL_DUPLETS_PERCENT["CG"] = 0
AL_DUPLETS_PERCENT["CA"] = 0
AL_DUPLETS_PERCENT["AC"] = 0
AL_DUPLETS_PERCENT["TG"] = 0
AL_DUPLETS_PERCENT["GC"] = 0



CA_TOTAL = 0
CA_DUPLETS = {}
CA_DUPLETS["TC"] = 0
CA_DUPLETS["AG"] = 0
CA_DUPLETS["GA"] = 0
CA_DUPLETS["CT"] = 0
CA_DUPLETS["AT"] = 0
CA_DUPLETS["GT"] = 0
CA_DUPLETS["TA"] = 0
CA_DUPLETS["CG"] = 0
CA_DUPLETS["CA"] = 0
CA_DUPLETS["AC"] = 0
CA_DUPLETS["TG"] = 0
CA_DUPLETS["GC"] = 0

CA_DUPLETS_PERCENT = {}
CA_DUPLETS_PERCENT["TC"] = 0
CA_DUPLETS_PERCENT["AG"] = 0
CA_DUPLETS_PERCENT["GA"] = 0
CA_DUPLETS_PERCENT["CT"] = 0
CA_DUPLETS_PERCENT["AT"] = 0
CA_DUPLETS_PERCENT["GT"] = 0
CA_DUPLETS_PERCENT["TA"] = 0
CA_DUPLETS_PERCENT["CG"] = 0
CA_DUPLETS_PERCENT["CA"] = 0
CA_DUPLETS_PERCENT["AC"] = 0
CA_DUPLETS_PERCENT["TG"] = 0
CA_DUPLETS_PERCENT["GC"] = 0



for sample in WM_SAMPLE_NAMES:

    SAMPLE_DUPLET_INPUT = f"{SPRINT_PROCESSED_OUTPUT_PATH}/{sample}/{sample}.ALL.SPRINT.duplets.tsv"

    duplets_file = open(SAMPLE_DUPLET_INPUT)
    read_duplets_file = csv.reader(duplets_file, delimiter="\t")

    dfile = []

    for row in read_duplets_file:
        dfile.append(row)

    body = dfile[0:]

    for row in range(0, len(body)):
        WM_DUPLETS[body[row][0]] += int(body[row][1])
        WM_TOTAL += int(body[row][1])
    
    
WM_DUPLETS_PERCENT["TC"] = format(round(((WM_DUPLETS["TC"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["AG"] = format(round(((WM_DUPLETS["AG"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["GA"] = format(round(((WM_DUPLETS["GA"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["CT"] = format(round(((WM_DUPLETS["CT"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["AT"] = format(round(((WM_DUPLETS["AT"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["GT"] = format(round(((WM_DUPLETS["GT"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["TA"] = format(round(((WM_DUPLETS["TA"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["CG"] = format(round(((WM_DUPLETS["CG"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["CA"] = format(round(((WM_DUPLETS["CA"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["AC"] = format(round(((WM_DUPLETS["AC"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["TG"] = format(round(((WM_DUPLETS["TG"] / WM_TOTAL)*100), 2), '.2f')
WM_DUPLETS_PERCENT["GC"] = format(round(((WM_DUPLETS["GC"] / WM_TOTAL)*100), 2), '.2f')    
    




for sample in NAWM_SAMPLE_NAMES:

    SAMPLE_DUPLET_INPUT = f"{SPRINT_PROCESSED_OUTPUT_PATH}/{sample}/{sample}.ALL.SPRINT.duplets.tsv"

    duplets_file = open(SAMPLE_DUPLET_INPUT)
    read_duplets_file = csv.reader(duplets_file, delimiter="\t")

    dfile = []

    for row in read_duplets_file:
        dfile.append(row)

    body = dfile[0:]

    for row in range(0, len(body)):
        NAWM_DUPLETS[body[row][0]] += int(body[row][1])
        NAWM_TOTAL += int(body[row][1])
    
    
NAWM_DUPLETS_PERCENT["TC"] = format(round(((NAWM_DUPLETS["TC"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["AG"] = format(round(((NAWM_DUPLETS["AG"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["GA"] = format(round(((NAWM_DUPLETS["GA"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["CT"] = format(round(((NAWM_DUPLETS["CT"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["AT"] = format(round(((NAWM_DUPLETS["AT"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["GT"] = format(round(((NAWM_DUPLETS["GT"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["TA"] = format(round(((NAWM_DUPLETS["TA"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["CG"] = format(round(((NAWM_DUPLETS["CG"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["CA"] = format(round(((NAWM_DUPLETS["CA"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["AC"] = format(round(((NAWM_DUPLETS["AC"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["TG"] = format(round(((NAWM_DUPLETS["TG"] / NAWM_TOTAL)*100), 2), '.2f')
NAWM_DUPLETS_PERCENT["GC"] = format(round(((NAWM_DUPLETS["GC"] / NAWM_TOTAL)*100), 2), '.2f')   




for sample in AL_SAMPLE_NAMES:

    SAMPLE_DUPLET_INPUT = f"{SPRINT_PROCESSED_OUTPUT_PATH}/{sample}/{sample}.ALL.SPRINT.duplets.tsv"

    duplets_file = open(SAMPLE_DUPLET_INPUT)
    read_duplets_file = csv.reader(duplets_file, delimiter="\t")

    dfile = []

    for row in read_duplets_file:
        dfile.append(row)

    body = dfile[0:]

    for row in range(0, len(body)):
        AL_DUPLETS[body[row][0]] += int(body[row][1])
        AL_TOTAL += int(body[row][1])
    
    
AL_DUPLETS_PERCENT["TC"] = format(round(((AL_DUPLETS["TC"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["AG"] = format(round(((AL_DUPLETS["AG"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["GA"] = format(round(((AL_DUPLETS["GA"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["CT"] = format(round(((AL_DUPLETS["CT"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["AT"] = format(round(((AL_DUPLETS["AT"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["GT"] = format(round(((AL_DUPLETS["GT"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["TA"] = format(round(((AL_DUPLETS["TA"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["CG"] = format(round(((AL_DUPLETS["CG"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["CA"] = format(round(((AL_DUPLETS["CA"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["AC"] = format(round(((AL_DUPLETS["AC"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["TG"] = format(round(((AL_DUPLETS["TG"] / AL_TOTAL)*100), 2), '.2f')
AL_DUPLETS_PERCENT["GC"] = format(round(((AL_DUPLETS["GC"] / AL_TOTAL)*100), 2), '.2f')   





for sample in CA_SAMPLE_NAMES:

    SAMPLE_DUPLET_INPUT = f"{SPRINT_PROCESSED_OUTPUT_PATH}/{sample}/{sample}.ALL.SPRINT.duplets.tsv"

    duplets_file = open(SAMPLE_DUPLET_INPUT)
    read_duplets_file = csv.reader(duplets_file, delimiter="\t")

    dfile = []

    for row in read_duplets_file:
        dfile.append(row)

    body = dfile[0:]

    for row in range(0, len(body)):
        CA_DUPLETS[body[row][0]] += int(body[row][1])
        CA_TOTAL += int(body[row][1])
    
    
CA_DUPLETS_PERCENT["TC"] = format(round(((CA_DUPLETS["TC"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["AG"] = format(round(((CA_DUPLETS["AG"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["GA"] = format(round(((CA_DUPLETS["GA"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["CT"] = format(round(((CA_DUPLETS["CT"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["AT"] = format(round(((CA_DUPLETS["AT"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["GT"] = format(round(((CA_DUPLETS["GT"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["TA"] = format(round(((CA_DUPLETS["TA"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["CG"] = format(round(((CA_DUPLETS["CG"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["CA"] = format(round(((CA_DUPLETS["CA"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["AC"] = format(round(((CA_DUPLETS["AC"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["TG"] = format(round(((CA_DUPLETS["TG"] / CA_TOTAL)*100), 2), '.2f')
CA_DUPLETS_PERCENT["GC"] = format(round(((CA_DUPLETS["GC"] / CA_TOTAL)*100), 2), '.2f')   



WM_duplets = ['T>C', 'A>G', 'G>A', 'C>T', 'A>T', 'G>T', 'T>A', 'C>G', 'C>A', 'A>C', 'T>G', 'G>C']
WM_duplets_values = list(WM_DUPLETS_PERCENT.values())
WM_duplets_values = [float(i) for i in WM_duplets_values]

NAWM_duplets = ['T>C', 'A>G', 'G>A', 'C>T', 'A>T', 'G>T', 'T>A', 'C>G', 'C>A', 'A>C', 'T>G', 'G>C']
NAWM_duplets_values = list(NAWM_DUPLETS_PERCENT.values())
NAWM_duplets_values = [float(i) for i in NAWM_duplets_values]

AL_duplets = ['T>C', 'A>G', 'G>A', 'C>T', 'A>T', 'G>T', 'T>A', 'C>G', 'C>A', 'A>C', 'T>G', 'G>C']
AL_duplets_values = list(AL_DUPLETS_PERCENT.values())
AL_duplets_values = [float(i) for i in AL_duplets_values]

CA_duplets = ['T>C', 'A>G', 'G>A', 'C>T', 'A>T', 'G>T', 'T>A', 'C>G', 'C>A', 'A>C', 'T>G', 'G>C']
CA_duplets_values = list(CA_DUPLETS_PERCENT.values())
CA_duplets_values = [float(i) for i in CA_duplets_values]


fig, axis = plt.subplots(2, 2, figsize=(15,15))

axis[0, 0].bar(WM_duplets, WM_duplets_values, color = "blue", bottom = 0)
axis[0, 0].set_title("CONTROL WM")
axis[0, 0].set_xlabel("Types of RNA editing sites", fontsize=12)
axis[0, 0].set_ylabel("Percent of RNA editing types (%)", fontsize=12)

axis[0, 1].bar(NAWM_duplets, NAWM_duplets_values, color = "darkblue", bottom = 0)
axis[0, 1].set_title("MS NAWM")
axis[0, 1].set_xlabel("Types of RNA editing sites", fontsize=12)
axis[0, 1].set_ylabel("Percent of RNA editing types (%)", fontsize=12)


axis[1, 0].bar(AL_duplets, AL_duplets_values, color = "darkviolet", bottom = 0)
axis[1, 0].set_title("MS AL")
axis[1, 0].set_xlabel("Types of RNA editing sites", fontsize=12)
axis[1, 0].set_ylabel("Percent of RNA editing types (%)", fontsize=12)


axis[1, 1].bar(CA_duplets, CA_duplets_values, color = "deepskyblue", bottom = 0)
axis[1, 1].set_title("MS CA")
axis[1, 1].set_xlabel("Types of RNA editing sites", fontsize=12)
axis[1, 1].set_ylabel("Percent of RNA editing types (%)", fontsize=12)

plt.suptitle("Distribution of 12 types of RNA editing events detected by SPRINT. Data represent the percentage of inferred RNA editing events detected in different sample types.", wrap=True, fontsize=16)

plt.savefig(f"{OUTPUT_POSTPROCESSING_PATH}/10_DUPLETS/MS_GSE138614_DUPLETS.pdf")

plt.show()






