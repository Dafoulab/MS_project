import csv
import subprocess

# mkdir RNAEditing_Database
# cd RNAEditing_Database
# wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
# gunzip TABLE1_hg38.txt.gz
# cut -f1,2,3,4,6,7,10 TABLE1_hg38.txt > REDIPortal_hg38_Aug_2022_Processed.tsv

INPUT_REDIPORTAL_PROCESSED_TSV_FILE = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS/RNAEditing_Database/REDIPortal_hg38_Aug_2022_Processed.tsv"
OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE = "/media/kokyriakidis/PRODUCTION/EDITING/SPRINT/MS/RNAEditing_Database/REDIPortal_hg38_Aug_2022_Processed.vcf"

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

with open(f"{OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}", mode='w') as vcfcandidates_file:
        vcfcandidates_file = csv.writer(vcfcandidates_file, delimiter='\t', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        vcfcandidates_file.writerow(['##fileformat=VCFv4.2'])
        vcfcandidates_file.writerow(['##reference=GRCh38/hg38'])
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
        vcfcandidates_file.writerow(['##INFO=<ID=RNAEDIT,Number=3,Type=String,Description="A known or predicted RNA-editing site",Source="REDIPortal",Version="Aug-2022">'])
        vcfcandidates_file.writerow(['##INFO=<ID=EDTYPE,Number=3,Type=String,Description="Type of a known or predicted RNA-editing site",Source="REDIPortal",Version="Aug-2022">'])
        vcfcandidates_file.writerow(['##INFO=<ID=EDREGION,Number=3,Type=String,Description="Functional region of a known or predicted RNA-editing site",Source="REDIPortal",Version="Aug-2022">'])
        vcfcandidates_file.writerow(['#CHROM', 'POS', 'ID',  'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
        with open(f"{INPUT_REDIPORTAL_PROCESSED_TSV_FILE}") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if row[0] == "Region":
                    next
                elif row[0] not in hg38_chr_names_to_include:
                    next
                else:
                    infofield = 'RNAEDIT=' + row[4] + ';EDTYPE=' + row[5] + ';EDREGION=' + row[6]
                    vcfcandidates_file.writerow([row[0], row[1], '.', row[2], row[3], '.', '.', infofield])


sorted_output = f"{OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}"[:-4] + '-sorted.vcf'
cmd = f"bcftools sort {OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE} > {sorted_output}"
subprocess.run(cmd, shell=True)

cmd = f"rm {OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}"
subprocess.run(cmd, shell=True)

cmd = f"rm {INPUT_REDIPORTAL_PROCESSED_TSV_FILE}"
subprocess.run(cmd, shell=True)

cmd = f"mv {sorted_output} {OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}"
subprocess.run(cmd, shell=True)

cmd = f"bgzip {OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}"
subprocess.run(cmd, shell=True)

cmd = f"bcftools index {OUTPUT_REDIPORTAL_PROCESSED_VCF_FILE}" + ".gz"
subprocess.run(cmd, shell=True)

