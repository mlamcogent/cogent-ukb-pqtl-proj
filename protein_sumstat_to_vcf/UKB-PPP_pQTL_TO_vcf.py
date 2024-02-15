#!/usr/bin/env python3

''''
This Python script is designed to convert the UKB-PPP Pqtl summary statistics file into the VCF (Variant Call Format) format. using gwas2vcf ; https://github.com/MRCIEU/gwas2vcf
'''
import pandas as pd
import numpy as np
import json,os
import argparse

parser=argparse.ArgumentParser(description="This Python script is designed to convert the UKB-PPP Pqtl summary statistics file into the VCF format.")
parser.add_argument('-filename','--filename', help="file name ", required=True)

args=parser.parse_args()
file=args.filename

path=os.getcwd()

df=pd.read_csv(file,sep="\t")
df=df[df['LOG10P']>=5]
df['LOG10P']=np.where(df['LOG10P']<200,df['LOG10P'],199)
df=df[['CHROM', 'GENPOS','ID','ALLELE1','ALLELE0','BETA','SE','N' ,'LOG10P','A1FREQ','INFO','GeneSymbol']]
df['LOG10P']=np.power(10,-df['LOG10P'])
ID=df.iloc[0,11]
df['ID']=df['ID'].str.replace(":","_")
df['CHROM']=np.where(df['CHROM']==23,"X",df['CHROM'])
df.to_csv(f'{file.split("/")[-1][:-3]}',sep="\t",index=None)

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
    "imp_info_col": 10, 
    "delimiter": "\t", 
    "header": "true", 
    "build": "GRCh38"}

with open('UKB-PPP_params.txt', 'w') as f:
  json.dump(paramsdict, f)

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/main.py \
    --out {path}/{file.split("/")[-1][:-7]}.vcf \
    --data {path}/{file.split("/")[-1][:-3]} \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz \
    --json UKB-PPP_params.txt \
    --id {ID} ''' )
os.system(f'tabix -f -p vcf  {path}/{file.split("/")[-1][:-7]}.vcf.gz')
