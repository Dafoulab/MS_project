##############################
##############################
### MS GSM4114296 ANALYSIS ###
##############################
##############################

import subprocess
import seaborn
import pandas as pd
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

OUTPUT_POSTPROCESSING_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/postprocessing"


##################################
#### Create directory structure ##
##################################
cmd = f"mkdir -p {OUTPUT_POSTPROCESSING_PATH}/11_EDITING_LEVELS"
subprocess.run(cmd, shell=True)



dataframe = pd.read_csv(f"{OUTPUT_POSTPROCESSING_PATH}/7_T_TEST/MS_GSM4114296_WMvsNAWMvsALvsCA_A2IC2U_SPRINT_QC_MICROGLIA_SAMPLE_Filtered_Plus_MeanAF_T_Test.tsv", sep='\t', header=0, low_memory=False)

d = {'WM_MEAN_FREQUENCY': dataframe.WM_MEAN_FREQUENCY.apply(lambda x: x*100), 'NAWM_MEAN_FREQUENCY':dataframe.NAWM_MEAN_FREQUENCY.apply(lambda x: x*100), 'AL_MEAN_FREQUENCY':dataframe.AL_MEAN_FREQUENCY.apply(lambda x: x*100), 'CA_MEAN_FREQUENCY':dataframe.CA_MEAN_FREQUENCY.apply(lambda x: x*100)}
df = pd.DataFrame(d)
seaborn.set(style = 'whitegrid')    
order = ['WM_MEAN_FREQUENCY', 'NAWM_MEAN_FREQUENCY', 'AL_MEAN_FREQUENCY', 'CA_MEAN_FREQUENCY']
ax = seaborn.violinplot(data=df[['WM_MEAN_FREQUENCY', 'NAWM_MEAN_FREQUENCY', 'AL_MEAN_FREQUENCY', 'CA_MEAN_FREQUENCY']], cut=0)
ax.set_xticklabels(['WM', 'NAWM', 'AL', 'CA'], fontsize=12)
ax.grid(False)
plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
plt.ylabel("RNA editing level (%)", size=12)


pairs=[("WM_MEAN_FREQUENCY", "NAWM_MEAN_FREQUENCY"), ("WM_MEAN_FREQUENCY", "AL_MEAN_FREQUENCY"), ("WM_MEAN_FREQUENCY", "CA_MEAN_FREQUENCY"), ("NAWM_MEAN_FREQUENCY", "AL_MEAN_FREQUENCY"), ("NAWM_MEAN_FREQUENCY", "CA_MEAN_FREQUENCY"), ("AL_MEAN_FREQUENCY", "CA_MEAN_FREQUENCY")]
annotator = Annotator(ax, pairs, data=df[['WM_MEAN_FREQUENCY', 'NAWM_MEAN_FREQUENCY', 'AL_MEAN_FREQUENCY', 'CA_MEAN_FREQUENCY']], order=order)
annotator.configure(test='Wilcoxon', comparisons_correction="Bonferroni", text_format='star', loc='inside')
annotator.apply_and_annotate()

plt.title("The distribution of RNA editing levels of MS sample types. \nStatistical significance was determined by Wilcoxon's signed-rank test and adjusted using the Bonferroni method. ****P < 0.001", wrap=True)

plt.savefig(f'{OUTPUT_POSTPROCESSING_PATH}/11_EDITING_LEVELS/MS_GSM4114296_EDITING_LEVELS.png')

plt.show()











