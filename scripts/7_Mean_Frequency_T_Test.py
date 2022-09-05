#############################
#############################
### MS GSE138614 ANALYSIS ###
#############################
#############################

import csv
import subprocess
from scipy.stats import ttest_ind

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSE138614/postprocessing"


##################################
### Create directory structure ###
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/7_T_TEST"
subprocess.run(cmd, shell=True)


##################################################################
### T TEST USING FILE AFTER QC, MICROGLIA AND SAMPLE FILTERING ###
##################################################################

variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/6_SAMPLE_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()


CA_MEAN_FREQUENCY_INDEX = header.index("CA_MEAN_FREQUENCY")

WM_FREQUENCY_SAMPLES = ["C3_1_FREQUENCY", "C5_8_FREQUENCY", "C5_10_FREQUENCY", "C5_13_FREQUENCY", "C4_16_FREQUENCY", "C2_18_FREQUENCY", "C4_19_FREQUENCY", "C1_24_FREQUENCY", "C5_25_FREQUENCY", "C3_26_FREQUENCY", "C4_28_FREQUENCY", "C2_32_FREQUENCY", "C5_33_FREQUENCY", "C4_34_FREQUENCY", "C3_53_FREQUENCY", "C4_54_FREQUENCY", "C3_55_FREQUENCY", "C1_56_FREQUENCY", "C1_57_FREQUENCY", "C2_58_FREQUENCY", "C1_59_FREQUENCY", "C1_60_FREQUENCY", "C2_61_FREQUENCY", "C2_62_FREQUENCY", "C3_63_FREQUENCY"]
NAWM_FREQUENCY_SAMPLES = ["MS3_11_FREQUENCY", "MS5_15_FREQUENCY", "MS8_20_FREQUENCY", "MS1_23_FREQUENCY", "MS3_27_FREQUENCY", "MS8_30_FREQUENCY", "MS3_31_FREQUENCY", "MS5_35_FREQUENCY", "MS3_41_FREQUENCY", "MS4_42_FREQUENCY", "MS1_43_FREQUENCY", "MS3_44_FREQUENCY", "MS1_45_FREQUENCY", "MS8_46_FREQUENCY", "MS6_47_FREQUENCY", "MS9_48_FREQUENCY", "MS8_49_FREQUENCY", "MS5_50_FREQUENCY", "MS6_51_FREQUENCY", "MS4_52_FREQUENCY", "MS5_94_FREQUENCY"]
AL_FREQUENCY_SAMPLES = ["MS5_64_FREQUENCY", "MS4_65_FREQUENCY", "MS6_66_FREQUENCY", "MS6_67_FREQUENCY", "MS5_68_FREQUENCY", "MS1_71_FREQUENCY", "MS6_72_FREQUENCY", "MS10_73_FREQUENCY", "MS1_74_FREQUENCY", "MS1_76_FREQUENCY", "MS2_77_FREQUENCY", "MS10_79_FREQUENCY", "MS7_86_FREQUENCY", "MS6_89_FREQUENCY", "MS4_90_FREQUENCY", "MS7_91_FREQUENCY"]
CA_FREQUENCY_SAMPLES = ["MS7_2_FREQUENCY", "MS6_3_FREQUENCY", "MS4_4_FREQUENCY", "MS8_5_FREQUENCY", "MS8_6_FREQUENCY", "MS8_9_FREQUENCY", "MS3_12_FREQUENCY", "MS4_14_FREQUENCY", "MS3_17_FREQUENCY", "MS4_21_FREQUENCY", "MS4_22_FREQUENCY", "MS4_29_FREQUENCY", "MS3_36_FREQUENCY", "MS8_37_FREQUENCY", "MS6_38_FREQUENCY", "MS10_39_FREQUENCY", "MS8_40_FREQUENCY"]

WM_FREQUENCY_INDICES = []
NAWM_FREQUENCY_INDICES = []
AL_FREQUENCY_INDICES = []
CA_FREQUENCY_INDICES = []

for cell in range(0, len(header)):
    if header[cell] in WM_FREQUENCY_SAMPLES:
        WM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in NAWM_FREQUENCY_SAMPLES:
        NAWM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in AL_FREQUENCY_SAMPLES:
        AL_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in CA_FREQUENCY_SAMPLES:
        CA_FREQUENCY_INDICES.append(int(cell))
    

for row in range(0, len(body)):
    WM_FREQUENCY_VALUES = []
    NAWM_FREQUENCY_VALUES = []
    AL_FREQUENCY_VALUES = []
    CA_FREQUENCY_VALUES = []
    for WM_FREQUENCY_INDEX in WM_FREQUENCY_INDICES:
        WM_FREQUENCY_VALUES.append(float(body[row][WM_FREQUENCY_INDEX]))
    for NAWM_FREQUENCY_INDEX in NAWM_FREQUENCY_INDICES:
        NAWM_FREQUENCY_VALUES.append(float(body[row][NAWM_FREQUENCY_INDEX]))
    for AL_FREQUENCY_INDEX in AL_FREQUENCY_INDICES:
        AL_FREQUENCY_VALUES.append(float(body[row][AL_FREQUENCY_INDEX]))
    for CA_FREQUENCY_INDEX in CA_FREQUENCY_INDICES:
        CA_FREQUENCY_VALUES.append(float(body[row][CA_FREQUENCY_INDEX]))
    WMvsALstat, WMvsALpval = ttest_ind(WM_FREQUENCY_VALUES, AL_FREQUENCY_VALUES)
    WMvsCAstat, WMvsCApval = ttest_ind(WM_FREQUENCY_VALUES, CA_FREQUENCY_VALUES)
    NAWMvsALstat, NAWMvsALpval = ttest_ind(NAWM_FREQUENCY_VALUES, AL_FREQUENCY_VALUES)
    NAWMvsCAstat, NAWMvsCApval = ttest_ind(NAWM_FREQUENCY_VALUES, CA_FREQUENCY_VALUES)
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+1), str(round(WMvsALpval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+2), str(round(WMvsCApval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+3), str(round(NAWMvsALpval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+4), str(round(NAWMvsCApval, 4)) )

header.insert((CA_MEAN_FREQUENCY_INDEX+1), "WMvsAL_T_TEST_PVALUE") 
header.insert((CA_MEAN_FREQUENCY_INDEX+2), "WMvsCA_T_TEST_PVALUE")
header.insert((CA_MEAN_FREQUENCY_INDEX+3), "NAWMvsAL_T_TEST_PVALUE")
header.insert((CA_MEAN_FREQUENCY_INDEX+4), "NAWMvsCA_T_TEST_PVALUE")  

fvfile = []
fvfile.append(header)
for row in range(0, len(body)):
    fvfile.append(body[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()






##########################################################
### T TEST USING FILE AFTER QC AND MICROGLIA FILTERING ###
##########################################################


variant_file = open(f"{OUTPUT_POSTPROCESSING_PATH}/5_MICROGLIA_FILTRATION/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_Filtered.tsv")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

variant_file.close()


CA_MEAN_FREQUENCY_INDEX = header.index("CA_MEAN_FREQUENCY")

WM_FREQUENCY_SAMPLES = ["C3_1_FREQUENCY", "C5_8_FREQUENCY", "C5_10_FREQUENCY", "C5_13_FREQUENCY", "C4_16_FREQUENCY", "C2_18_FREQUENCY", "C4_19_FREQUENCY", "C1_24_FREQUENCY", "C5_25_FREQUENCY", "C3_26_FREQUENCY", "C4_28_FREQUENCY", "C2_32_FREQUENCY", "C5_33_FREQUENCY", "C4_34_FREQUENCY", "C3_53_FREQUENCY", "C4_54_FREQUENCY", "C3_55_FREQUENCY", "C1_56_FREQUENCY", "C1_57_FREQUENCY", "C2_58_FREQUENCY", "C1_59_FREQUENCY", "C1_60_FREQUENCY", "C2_61_FREQUENCY", "C2_62_FREQUENCY", "C3_63_FREQUENCY"]
NAWM_FREQUENCY_SAMPLES = ["MS3_11_FREQUENCY", "MS5_15_FREQUENCY", "MS8_20_FREQUENCY", "MS1_23_FREQUENCY", "MS3_27_FREQUENCY", "MS8_30_FREQUENCY", "MS3_31_FREQUENCY", "MS5_35_FREQUENCY", "MS3_41_FREQUENCY", "MS4_42_FREQUENCY", "MS1_43_FREQUENCY", "MS3_44_FREQUENCY", "MS1_45_FREQUENCY", "MS8_46_FREQUENCY", "MS6_47_FREQUENCY", "MS9_48_FREQUENCY", "MS8_49_FREQUENCY", "MS5_50_FREQUENCY", "MS6_51_FREQUENCY", "MS4_52_FREQUENCY", "MS5_94_FREQUENCY"]
AL_FREQUENCY_SAMPLES = ["MS5_64_FREQUENCY", "MS4_65_FREQUENCY", "MS6_66_FREQUENCY", "MS6_67_FREQUENCY", "MS5_68_FREQUENCY", "MS1_71_FREQUENCY", "MS6_72_FREQUENCY", "MS10_73_FREQUENCY", "MS1_74_FREQUENCY", "MS1_76_FREQUENCY", "MS2_77_FREQUENCY", "MS10_79_FREQUENCY", "MS7_86_FREQUENCY", "MS6_89_FREQUENCY", "MS4_90_FREQUENCY", "MS7_91_FREQUENCY"]
CA_FREQUENCY_SAMPLES = ["MS7_2_FREQUENCY", "MS6_3_FREQUENCY", "MS4_4_FREQUENCY", "MS8_5_FREQUENCY", "MS8_6_FREQUENCY", "MS8_9_FREQUENCY", "MS3_12_FREQUENCY", "MS4_14_FREQUENCY", "MS3_17_FREQUENCY", "MS4_21_FREQUENCY", "MS4_22_FREQUENCY", "MS4_29_FREQUENCY", "MS3_36_FREQUENCY", "MS8_37_FREQUENCY", "MS6_38_FREQUENCY", "MS10_39_FREQUENCY", "MS8_40_FREQUENCY"]

WM_FREQUENCY_INDICES = []
NAWM_FREQUENCY_INDICES = []
AL_FREQUENCY_INDICES = []
CA_FREQUENCY_INDICES = []

for cell in range(0, len(header)):
    if header[cell] in WM_FREQUENCY_SAMPLES:
        WM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in NAWM_FREQUENCY_SAMPLES:
        NAWM_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in AL_FREQUENCY_SAMPLES:
        AL_FREQUENCY_INDICES.append(int(cell))
    elif header[cell] in CA_FREQUENCY_SAMPLES:
        CA_FREQUENCY_INDICES.append(int(cell))
    

for row in range(0, len(body)):
    WM_FREQUENCY_VALUES = []
    NAWM_FREQUENCY_VALUES = []
    AL_FREQUENCY_VALUES = []
    CA_FREQUENCY_VALUES = []
    for WM_FREQUENCY_INDEX in WM_FREQUENCY_INDICES:
        WM_FREQUENCY_VALUES.append(float(body[row][WM_FREQUENCY_INDEX]))
    for NAWM_FREQUENCY_INDEX in NAWM_FREQUENCY_INDICES:
        NAWM_FREQUENCY_VALUES.append(float(body[row][NAWM_FREQUENCY_INDEX]))
    for AL_FREQUENCY_INDEX in AL_FREQUENCY_INDICES:
        AL_FREQUENCY_VALUES.append(float(body[row][AL_FREQUENCY_INDEX]))
    for CA_FREQUENCY_INDEX in CA_FREQUENCY_INDICES:
        CA_FREQUENCY_VALUES.append(float(body[row][CA_FREQUENCY_INDEX]))
    WMvsALstat, WMvsALpval = ttest_ind(WM_FREQUENCY_VALUES, AL_FREQUENCY_VALUES)
    WMvsCAstat, WMvsCApval = ttest_ind(WM_FREQUENCY_VALUES, CA_FREQUENCY_VALUES)
    NAWMvsALstat, NAWMvsALpval = ttest_ind(NAWM_FREQUENCY_VALUES, AL_FREQUENCY_VALUES)
    NAWMvsCAstat, NAWMvsCApval = ttest_ind(NAWM_FREQUENCY_VALUES, CA_FREQUENCY_VALUES)
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+1), str(round(WMvsALpval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+2), str(round(WMvsCApval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+3), str(round(NAWMvsALpval, 4)) )
    body[row].insert((CA_MEAN_FREQUENCY_INDEX+4), str(round(NAWMvsCApval, 4)) )

header.insert((CA_MEAN_FREQUENCY_INDEX+1), "WMvsAL_T_TEST_PVALUE") 
header.insert((CA_MEAN_FREQUENCY_INDEX+2), "WMvsCA_T_TEST_PVALUE")
header.insert((CA_MEAN_FREQUENCY_INDEX+3), "NAWMvsAL_T_TEST_PVALUE")
header.insert((CA_MEAN_FREQUENCY_INDEX+4), "NAWMvsCA_T_TEST_PVALUE")  

fvfile = []
fvfile.append(header)
for row in range(0, len(body)):
    fvfile.append(body[row])

f = open(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSE138614_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_Filtered_Plus_MeanAF_T_Test.tsv", "w")

for line in fvfile:
    print("\t".join(line), file=f)

f.close()







































