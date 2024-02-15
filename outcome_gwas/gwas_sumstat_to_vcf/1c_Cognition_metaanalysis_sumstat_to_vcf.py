#!/usr/bin/env python3

''''
This Python script is designed to convert the meta-analyzed GWAS cognition summary statistics file (Cognition_Savage_iq_metaanalysis_noUKB and UKBB without UKB-PPP meta-analysis) to VCF format.
'''

import pandas as pd
import numpy as np
import json,os

filename="Cognition_Meta_GWAS_without_UKBPP1.tbl" #location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/Cognition/"
fdf=pd.read_csv(filename,sep="\t")
fdf=fdf[['MarkerName', 'Allele1', 'Allele2', 'Weight','Zscore', 'P-value','Freq1']]

#Calculate the beta and standard error (SE) values from the Z score  : Reference https://ctg.cncr.nl/documents/p1651/readme.txt ;  https://www.biostars.org/p/319584/
fdf['beta'] = fdf['Zscore'] / np.sqrt(2 * fdf['Freq1'] * (1 - fdf['Freq1']) * (fdf['Weight'] + fdf['Zscore']**2))
fdf['se']   = 1 / np.sqrt(2 * fdf['Freq1'] * (1 - fdf['Freq1']) * (fdf['Weight'] + fdf['Zscore']**2))

fdf[['CHR','BP']]=fdf["MarkerName"].str.split("_",expand=True)[[0,1]]
fdf=fdf.rename(columns={'MarkerName':'SNP','Weight':'totalN','P-value':'P'})
fdf=fdf[['CHR','SNP','BP', 'Allele1','Allele2','Freq1','beta','se','P','totalN']]

fdf['Allele1']=fdf['Allele1'].str.upper()
fdf['Allele2']=fdf['Allele2'].str.upper()

fdf3=fdf[~fdf["CHR"].isna()]
fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")
fdf3[['BP','totalN']]=fdf3[['BP','totalN']].astype("int")

fdf3[['Freq1', 'beta','se','P']]=fdf3[['Freq1', 'beta','se','P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','Allele1','Allele2','beta','se','totalN','P','Freq1']]
fdf3.to_csv(f'{filename[:-4]}.tsv',sep="\t",index=None)

paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "pval_col": 8,
    "eaf_col": 9,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh38"
    }

with open('Cognition_Meta_GWAS_without_UKBPP1_params.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format 
path=os.getcwd()
ID="Cognition_Meta"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --json Cognition_Meta_GWAS_without_UKBPP1_params.txt \
    --id {ID} > {filename[:-4]}.error 2>&1 ''' )


os.system(f"mv {path}/{filename[:-4]}.vcf.gz {filename[:-4]}_GRCh38.vcf.gz")
os.system(f"mv {path}/{filename[:-4]}.vcf.gz.tbi {filename[:-4]}_GRCh38.vcf.gz.tbi")

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-4]}_GRCh38.vcf.gz | bgzip -c > {path}/{filename[:-4]}_GRCh38_sort.vcf.gz")
os.system(f'tabix -f -p vcf {path}/{filename[:-4]}_GRCh38_sort.vcf.gz')



##Create unique ID using CHR, POS, REF AND ALT columns
os.system(f'zgrep -v "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38.tab')
os.system(f'zgrep  "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38_header.txt')

df=pd.read_csv(f"{path}/{filename[:-4]}_GRCh38.tab",sep="\t")
format_df=df[ID].str.split(":",expand=True)
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]
df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')

os.system(f" rm {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab {path}/{filename[:-4]}_GRCh38.tab  ")
