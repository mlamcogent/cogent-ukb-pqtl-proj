#!/usr/bin/env python3

''''
This Python script is designed to convert the daner_bip_pgc3_nm_noukbiobank.gz summary statistics file to VCF format.
'''

import pandas as pd
import numpy as np
import json,os


##Downloaded from PGC:  https://pgc.unc.edu/for-researchers/download-results/
os.system("wget https://figshare.com/ndownloader/files/40036705")
os.system("mv 40036705 daner_bip_pgc3_nm_noukbiobank.gz")
os.system('wget https://figshare.com/ndownloader/files/40036729')
os.system("mv 40036729 readme.txt")


filename="daner_bip_pgc3_nm_noukbiobank.gz" 

fdf=pd.read_csv(filename,sep="\t")
fdf2=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_313436', 'INFO','OR', 'SE', 'P','Nca', 'Nco']]
fdf3=fdf2[~fdf2["SNP"].isna()]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['FRQ_U_313436', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_313436', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3["Beta"]=np.log(fdf3["OR"])
fdf3.drop("OR",axis=1,inplace=True)
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"Beta",'SE','Nco','Nca','P','FRQ_U_313436','INFO' ]]
fdf3.to_csv(f'{filename[:-3]}.tsv',sep="\t",index=None)

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
    "build": "GRCh37"
    }

with open('daner_bip_pgc3_params.txt', 'w') as f:
  json.dump(paramsdict, f)


path=os.getcwd()

ID="BIP_PGC3_noukb"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-3]}.vcf \
    --data {path}/{filename[:-3]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json daner_bip_pgc3_params.txt \
    --id {ID} > {filename[:-3]}.error 2>&1 ''' )

os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}.vcf.gz')

##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename[:-3]}.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-3]}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-3]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-3]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38.vcf.gz')
os.system("rm daner_bip_pgc3_nm_noukbiobank_GRCh38.vcf ")

##Create unique ID using CHR, POS, REF AND ALT columns
os.system(f'zgrep -v "##" {path}/{filename[:-3]}_GRCh38.vcf.gz > {path}/{filename[:-3]}_GRCh38.tab')
os.system(f'zgrep  "##" {path}/{filename[:-3]}_GRCh38.vcf.gz > {path}/{filename[:-3]}_GRCh38_header.txt')

df=pd.read_csv(f"{path}/{filename[:-3]}_GRCh38.tab",sep="\t")
format_df=df["BIP_PGC3_noukb"].str.split(":",expand=True)
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[7]=df["ID"]
format_df=format_df.astype("str")
df['BIP_PGC3_noukb']=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df['FORMAT'].str.split(":",expand=True).isna().sum()
df['BIP_PGC3_noukb'].str.split(":",expand=True).isna().sum()
print(df)

df.to_csv(f"{path}/{filename[:-3]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-3]}_GRCh38_header.txt {path}/{filename[:-3]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-3]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38_UniqID.vcf.gz')
os.system(f" rm {path}/{filename[:-3]}_GRCh38_header.txt {path}/{filename[:-3]}_GRCh38_UniqID.tab {path}/{filename[:-3]}_GRCh38.tab  ")
os.system(f"rm {path}/{filename[:-3]}.vcf.gz* {path}/{filename[:-3]}_GRCh38.vcf.gz*")