!/bin/bash

UNTRIMMED_INPUT_FASTQ_FOLDER_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/SPRINT/input"
TRIMMED_FASTQ_FOLDER_PATH = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS_GSM4114296/SPRINT/input-trimmed"

cd ${UNTRIMMED_INPUT_FASTQ_FOLDER_PATH}

for f1 in *_1.fastq
do
       f2=${f1%%_1.fastq}"_2.fastq"
       f1_trimmed = ${f1%%_1.fastq}"-trimmed_1.fastq"
       f2_trimmed = ${f1%%_1.fastq}"-trimmed_2.fastq"
       fastp -i $f1 -I $f2 -q 25 -u 10 -l 50 -y -x -w 15 -o "${TRIMMED_FASTQ_FOLDER_PATH}/$f1_trimmed" -O "${TRIMMED_FASTQ_FOLDER_PATH}/$f2_trimmed"
done