#!/usr/bin/env python3


''''
This python script is to convert the "Cognition_Savage_iq_metaanalysis_noUKB.out.gz" sumstat file to vcf format, this will be later used for meta analysis

'''

import pandas as pd
import numpy as np
import json,os

sumstat_file="Cognition_Savage_iq_metaanalysis_noUKB.out.gz"
os.system('zgrep -v "##" {sumstat_file} > Cognition_Savage_iq_metaanalysis_noUKB.tsv')
filename="Cognition_Savage_iq_metaanalysis_noUKB.tsv"
fdf=pd.read_csv(filename,sep="\t")
fdf=fdf[['CHR','SNP','BP', 'Allele1', 'Allele2','totalN','Zscore', 'P', 'EAF_HRC']]

#Calculate the beta and standard error (SE) values from the Z score  : Reference https://ctg.cncr.nl/documents/p1651/readme.txt ;  https://www.biostars.org/p/319584/
fdf['beta'] = fdf['Zscore'] / np.sqrt(2 * fdf['EAF_HRC'] * (1 - fdf['EAF_HRC']) * (fdf['totalN'] + fdf['Zscore']**2))
fdf['se']   = 1 / np.sqrt(2 * fdf['EAF_HRC'] * (1 - fdf['EAF_HRC']) * (fdf['totalN'] + fdf['Zscore']**2))

## Select the relevant columns
fdf=fdf[['CHR','SNP','BP', 'Allele1','Allele2','EAF_HRC','beta','se','P','totalN']]
fdf['Allele1']=fdf['Allele1'].str.upper()
fdf['Allele2']=fdf['Allele2'].str.upper()

fdf3=fdf[~fdf["CHR"].isna()]
fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")
fdf3[['BP','totalN']]=fdf3[['BP','totalN']].astype("int")

fdf3[['EAF_HRC', 'beta','se','P']]=fdf3[['EAF_HRC', 'beta','se','P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','Allele1','Allele2','beta','se','totalN','P','EAF_HRC']]
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
    "build": "GRCh37"
    }

with open('Cognition_Savage_iq_metaanalysis_noUKB_params.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format using gwas2vcf : Reff https://github.com/MRCIEU/gwas2vcf
path=os.getcwd()
ID="Cognition_Savage_NoUKB"

os.system(f'''python Software/gwas2vcf/main.py \
                --out {path}/{filename[:-4]}.vcf \
                --data {path}/{filename[:-4]}.tsv \
                --ref gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                --json Cognition_Savage_iq_metaanalysis_noUKB_params.txt \
                --id {ID} > {filename[:-4]}.error 2>&1 ''' )



##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename[:-4]}.vcf.gz \
              gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-4]}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-4]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-4]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38.vcf.gz')
os.system(f"rm {path}/{filename[:-4]}_GRCh38.vcf ")


##Create unique ID using CHR, POS, REF AND ALT columns
os.system(f'zgrep -v "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38.tab')
os.system(f'zgrep  "##" {path}/{filename[:-4]}_GRCh38.vcf.gz > {path}/{filename[:-4]}_GRCh38_header.txt')

df=pd.read_csv(f"{path}/{filename[:-4]}_GRCh38.tab",sep="\t")


##Ading UniqID to the format section
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df=df[ID].str.split(":",expand=True)
format_df=format_df.astype("str")
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]
df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')


## For cgnition meta analysis the file require in tab separate format
df2=pd.concat([df.iloc[:,0:5],format_df],axis=1).iloc[:,:-1]

#The p-value is in -log10 format; for meta-analysis, convert this value to the normal p-value.
df2['LP']=np.power(10,-df2['LP'])
df2.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID_ForMeta.tsv",sep="\t",index=None)
