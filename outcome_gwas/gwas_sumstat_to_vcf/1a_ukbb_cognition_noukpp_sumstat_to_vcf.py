#!/usr/bin/env python3

''''
This python script is to convert the "UKB Cognitive GWAS without UKB-PPP samples" sumstat file to vcf format, this will be later used for meta analysis

'''

import pandas as pd
import numpy as np
import tempfile
import json,os


os.system('zgrep -v "##" regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt.gz > regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.tsv')
filename="regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.tsv"

fdf=pd.read_csv(filename,sep=" ")
fdf=fdf[['CHROM','ID','GENPOS', 'ALLELE0', 'ALLELE1','N', 'A1FREQ', 'INFO', 'BETA', 'SE','LOG10P']]

##Convert the -log10 P value to P valuee
fdf['LOG10P']=np.power(10,-fdf['LOG10P'])
fdf=fdf[['CHROM','ID','GENPOS','ALLELE0','ALLELE1','A1FREQ','BETA','SE','LOG10P','N','INFO']]

fdf3=fdf[~fdf["CHROM"].isna()]
fdf3["CHROM"]=fdf3["CHROM"].astype("int").astype("str")
fdf3[['GENPOS','N']]=fdf3[['GENPOS','N']].astype("int")
fdf3[['A1FREQ','BETA','SE','LOG10P']]=fdf3[['A1FREQ','BETA','SE','LOG10P']].astype("float")
fdf3=fdf3[['CHROM','GENPOS','ID','ALLELE1','ALLELE0','BETA','SE','N','LOG10P','A1FREQ','INFO']]

##Remove the variants with INFO score <0.6 & MAF 0.005
fdf3=fdf3[ ((fdf3['INFO']>=0.6 ) & (fdf3['A1FREQ']>=0.005 ) & (fdf3['A1FREQ']<0.995 ) )]
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
    "imp_info_col":10,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"
    }

with open('Cognition_UKB_params.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format 
path=os.getcwd()
ID="Cognition_UKB"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json Cognition_UKB_params.txt \
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
format_df[6]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')

## For cgnition meta analysis the file require in tab separate format
df2=pd.concat([df.iloc[:,0:5],format_df],axis=1).iloc[:,:-1]
df2.columns=['#CHROM','POS','ID','REF','ALT','ES', 'SE', 'LP', 'AF', 'SS', 'SI']
df2['LP']=df2['LP'].astype("float")
df2['LP']=np.power(10,-df2['LP'])
df2.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID_ForMeta.tsv",sep="\t",index=None)
os.system(f" rm {path}/{filename[:-4]}_GRCh38_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab {path}/{filename[:-4]}_GRCh38.tab  ")
