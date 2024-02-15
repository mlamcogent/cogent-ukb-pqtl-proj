#!/usr/bin/env python3

''''
This Python script is designed to convert the iPSYCH_Depression_excl_UKB_23andMe summary statistics file to VCF format.
'''

import pandas as pd
import numpy as np
import json,os



##File downloaded from https://ipsych.dk/en/research/downloads
os.system("wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz")
os.system("wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/Readme_excl_UKB.pdf")


filename="daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz"
os.system(f'''zgrep -v "##" {filename} >{filename[:-3]}.tsv''')
filename="daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.tsv" 
fdf=pd.read_csv(f"{filename[:-3]}.tsv",sep=" ")


fdf3=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_507679', 'INFO', 'OR', 'SE', 'P', 'Nca', 'Nco']]
fdf3=fdf3[~fdf3["CHR"].isna()]
fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")

fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_507679', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"OR",'SE','Nco','Nca','P','FRQ_U_507679','INFO' ]]

## Effect size column in OR format, so convert to beta format
fdf3["OR"]=np.log(fdf3["OR"])
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

with open('Depression_iPSYCH_2023_params.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format 
path=os.getcwd()
ID="Depression_iPSYCH_2023"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json Depression_iPSYCH_2023_params.txt \
    --id {ID} > {filename[:-4]}.error 2>&1 ''' )



##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename[:-4]}.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-4]}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-4]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-4]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38.vcf.gz')
os.system(f"rm {path}/{filename[:-4]}_GRCh38.vcf ")

##Create unique ID using CHR, POS, REF AND ALT columns
os.system(f'zgrep -v "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38.tab')
os.system(f'zgrep  "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38_header.txt')

df=pd.read_csv(f"{path}/{filename[:-4]}_GRCh38.tab",sep="\t")
format_df=df[ID].str.split(":",expand=True)
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[7]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')

os.system(f" rm {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab {path}/{filename[:-4]}_GRCh38.tab  ")