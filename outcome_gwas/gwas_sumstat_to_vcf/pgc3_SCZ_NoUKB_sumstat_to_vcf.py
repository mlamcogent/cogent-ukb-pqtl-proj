#!/usr/bin/env python3

''''
This Python script is designed to convert the daner_bip_pgc3_nm_noukbiobank.gz summary statistics file to VCF format.
'''

import pandas as pd
import numpy as np
import json,os


##This summary statistics file was obtained from PGC upon request
filename="daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios" 
fdf=pd.read_csv(filename,sep="\t")
fdf3=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_93456', 'INFO', 'OR', 'SE', 'P', 'Nca', 'Nco']]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")

#The effect size, in the file reported as odds ratio (OR), has been converted to beta.
fdf3["BETA"]=np.log(fdf3["OR"])
fdf3[['FRQ_U_93456', 'INFO', 'BETA', 'SE', 'P']]=fdf3[['FRQ_U_93456', 'INFO', 'BETA', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"BETA",'SE','Nco','Nca','P','FRQ_U_93456','INFO' ]]
fdf3.to_csv(f'{filename}.tsv',sep="\t",index=None)


paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "ncase_col":8,
    "pval_col": 9,
    "eaf_col": 10,
    "imp_info_col":11,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"}

with open('pgc3_scz_params.txt', 'w') as f:
  json.dump(paramsdict, f)


##Convert to vcf format 
path=os.getcwd()
ID="PGC3_SCZ_NoUKB"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename}.vcf \
    --data {path}/{filename}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json pgc3_scz_params.txt \
    --id {ID} > {filename}.error 2>&1 ''' )

#os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}.vcf.gz')

##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename}.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename}_GRCh38.vcf | bgzip -c > {path}/{filename}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename}_GRCh38.vcf.gz')
os.system(f"rm {path}/{filename}_GRCh38.vcf ")


##Create unique ID using CHR, POS, REF AND ALT columns
os.system(f'zgrep -v "##" {path}/{filename}_GRCh38.vcf.gz > {path}/{filename}_GRCh38.tab')
os.system(f'zgrep  "##" {path}/{filename}_GRCh38.vcf.gz > {path}/{filename}_GRCh38_header.txt')

df=pd.read_csv(f"{path}/{filename}_GRCh38.tab",sep="\t")
format_df=df[ID].str.split(":",expand=True)
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[7]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df.to_csv(f"{path}/{filename}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename}_GRCh38_header.txt {path}/{filename}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename}_GRCh38_UniqID.vcf.gz')
os.system(f" rm {path}/{filename}_GRCh38_header.txt {path}/{filename}_GRCh38_UniqID.tab {path}/{filename}_GRCh38.tab  ")