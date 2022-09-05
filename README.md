# RNA Editing Analysis of Multiple Sclerosis (MS) RNASeq samples 

This repository contains the code used to generate the figures and tables of the RNA editing analysis. 


## Environment setup

We recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to create a separate environment to avoid changing any installed software in your system. This is also the easiest way to install the required software. 

[SPRINT](https://github.com/jumphone/SPRINT) needs Python 2.7 installed as well as Samtools v1.2 and BWA v0.7.12. So, we will need to create 2 seperate conda environments, one to run SPRINT and one to do the postprocessing of SPRINT output. 

1. Install Miniconda

You can use the following command to install Miniconda. During installation, you may need to press enter and type "yes".

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
2. Create a conda environment to run SPRINT
   
Create a new conda environment named runSPRINT from the `runSPRINT.yml` file, found inside the `envs` folder, with the following command:

```bash
conda env create -f runSPRINT.yml
```

1. Create a conda environment to postprocess SPRINT output
   
Create a new conda environment named postprocessSPRINT from the `postprocessSPRINT.yml` file, found inside the `envs` folder, with the following command:

```bash
conda env create -f postprocessSPRINT.yml
```

## Bcbio setup

We used [bcbio v1.2.9](https://github.com/bcbio/bcbio-nextgen) to perform bulk RNASeq analysis and retrieve sequencing and mapping metrics. `bcbio` is a validated, scalable, community developed variant calling, RNA-seq and small RNA analysis analysis platform. 

The commands used to install `bcbio` can be found in the `bcbio_install_commands.sh` file inside the `bcbio/install` folder. Further installation instructions can be found [here](https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html).

The configuration file used to run the bcbio bulk RNASeq analysis, named `MS_GSM4114296_bcbio_config.yaml`, can be found inside the `bcbio/config` folder.

The `MS_GSM4114296_Sample_Metadata.csv` file, found inside the `bcbio/metadata` folder, contains the sample metadata.

The `MS_GSM4114296_Sample_Metrics.csv` file, found inside the `bcbio/metrics` folder, contains several sequencing and mapping metrics used to exclude samples from further RNA Editing analysis with SPRINT.


## Identification of RNA editing sites using SPRINT

#### Download Reference

We used the GRCh38 reference genome. You can download it using the following commands:

```bash
mkdir -p reference

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
```

#### FastQ Trimming

We used [fastp](https://github.com/OpenGene/fastp) to trim adapter sequences. Trimming commands can be found in the `MS_GSM4114296_Trimming_Commands.txt` file, inside the `SPRINT/scripts` folder. 

#### Run SPRINT

Instruction on how to download, prepare and run SPRINT can be found [here](https://github.com/jumphone/SPRINT).

The commands used to run SPRINT can be found in the `MS_GSM4114296_runSPRINT.sh` file, inside the `SPRINT/scripts` folder.

## Postprocess SPRINT output

#### Download REDIPortal V2.0 Database

We used the REDIPortal V2.0 RNA Editing Database. You can download it using the following commands:

```bash
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
gunzip TABLE1_hg38.txt.gz
cut -f1,2,3,4,6,7,10 TABLE1_hg38.txt > REDIPortal_hg38_Processed.tsv

```

Then, use the scipt `create_RNAEditing_db.py` inside the `scripts` folder to process it and convert it to a VCF format file.


#### Postprocessing Analysis

All scripts used to generate figures and tables are located inside the `scripts` folder.

All resources needed to run the scripts are located inside the `resources` folder.



















