
import os
import pandas as pd

''''
This Python script is designed to remove the samples used in pQTL analysis and then perform GWAS without those samples."
'''
os.system("module load PLINK/2.00a2.3_x86_64")
os.system("module load Regenie/3.2.4-GCC-11.2.0")
os.system("module load bgen/1.1.6")
os.system("module load BGEN-enkre/1.1.7-GCC-11.2.0")


##Samples with cognitive phenotype
cog_phenotype="/data/pheno_v3/ukbb500k.fiall_baseline_model2.EUR.regenie.v3.txt"
#Samples used in the pqtl analysis
UKB_PPP_fam="/data/ukb_imp_chr1_mac50_info07_b0_7_patched_bfiles.fam" 
# binary files of ukbiobank samples 
genotype_file="/data/ukb_chrauto_v2_snpqc"



## Due to the absence of family IDs and individual IDs in the binary file, the file needs to be converted to BED, BIM, and FAM formats to acquire this essential information.
os.system(f"plink2 --pfile {genotype_file} --make-bed --out ukb_chrauto_v2_snpqc")


#Identify the shared samples between the entire UK Biobank dataset (UKBB500K) and the samples used in the pQTL analysis (UKB_PPP).

UKB_PPP_fam_df=pd.read_csv(UKB_PPP_fam,sep="\t",header=None)
UKB_PPP_fam_df.columns=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"]

ukbb500k_fam_df=pd.read_csv("ukb_chrauto_v2_snpqc.fam",sep="\t",header=None)
ukbb500k_fam_df.columns=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex","ukbb_disease"]

##The samples present in both UKBB500K and UKB_PPP will not be included in the GWAS analysis
common_samples=pd.merge(ukbb500k_fam_df,UKB_PPP_fam_df,left_on=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex",
                            "ukbb_disease"], right_on=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"]).iloc[:,:2]
common_samples.to_csv("UKB_PPP_Samples_to_removefrom_CogGWAS.tsv",index=None,header=None)



##Remove samples from the phenotype file if they have been used in the UKB_PPP.
ukbb500k_UKB_PPP_df=pd.merge(ukbb500k_fam_df,UKB_PPP_fam_df,left_on=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex",
                            "ukbb_disease"], right_on=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"],how="outer")
ukbb500k_UKB_specific_df=ukbb500k_UKB_PPP_df[ukbb500k_UKB_PPP_df["ukpp_IID"].isna()].iloc[:,0:2]

##Read the cognitive phenotype file
cog_phenotype_df=pd.read_csv(cog_phenotype,sep=" ")
cog_phenotype_withoutUKB_PPP_df=pd.merge(ukbb500k_UKB_specific_df,cog_phenotype_df,left_on=["ukbb_FID","ukbb_IID"],
                                          right_on=["FID","IID"]).drop(["ukbb_FID","ukbb_IID"],axis=1)

##Save the cognitive phenotype file after removing the samples used in the UKB_PPP.
cog_phenotype_withoutUKB_PPP_df.to_csv("ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt",sep=" ",index=None)



# QC UKBB genotype data for Regenie based on the following code
os.system(f'''plink2 --threads 48 --pfile {genotype_file} \
            --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
            --remove UKB_PPP_Samples_to_removefrom_CogGWAS.tsv \
            --keep /data/ukbb500k.cogphenopc.EUR.feid.v1.keep \
            --write-snplist --write-samples --no-id-header --out ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR ''')
            

#-------------------
# Regenie Step 1: fitting the null genome-wide ridge regression
# model 2 : sex + age + PC + assessment center
#-------------------

os.system(f'''regenie \
        --step 1 \
        --pgen /data/Genotyped/ukb_chrauto_v2_snpqc \
        --extract ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.snplist \
        --keep ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.id \
        --covarFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
        --phenoFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
        --phenoColList fiall_baseline_irnt \
        --covarColList sex,fiall_baseline_age_mean0,fiall_baseline_age_mean02,fiall_baseline_age_mean0_sex,fiall_baseline_age_mean02_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,assess_center10003,assess_center11001,assess_center11002,assess_center11003,assess_center11004,assess_center11005,assess_center11006,assess_center11007,assess_center11008,assess_center11009,assess_center11011,assess_center11012,assess_center11013,assess_center11014,assess_center11016,assess_center11017,assess_center11018,assess_center11020,assess_center11021,assess_center11022,assess_center11023 \
        --bsize 1000 \
        --threads 48 \
        --cv 10 \
        --out regenie_ukb_step1.fiall_baseline_irnt_model2''')

#-------------------
# Regenie Step 2: performing single-variant association tests
#-------------------

for chr in range(1,23):
   os.system(f''' 
    regenie \
    --step 2 \
    --bgen /data/ukb_imp_chr{chr}_v3.bgen \
    --sample /data/ukb_imp_v3.sample \
    --keep ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.id \
    --covarFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
    --phenoFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
    --phenoColList fiall_baseline_irnt \
    --covarColList sex,fiall_baseline_age_mean0,fiall_baseline_age_mean02,fiall_baseline_age_mean0_sex,fiall_baseline_age_mean02_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,assess_center10003,assess_center11001,assess_center11002,assess_center11003,assess_center11004,assess_center11005,assess_center11006,assess_center11007,assess_center11008,assess_center11009,assess_center11011,assess_center11012,assess_center11013,assess_center11014,assess_center11016,assess_center11017,assess_center11018,assess_center11020,assess_center11021,assess_center11022,assess_center11023 \
    --bsize 1000  \
    --pred regenie_ukb_step1.fiall_baseline_irnt_model2_pred.list \
    --chr {chr} --threads 48 \
    --out regenie_ukb_step2_linear_model2_chr{chr}''')
 

#------------
# Merge chromosome-wise Regenie GWAS output into a single file.
#-------------------

os.system(f''' cp regenie_ukb_step2_linear_model2_chr1_fiall_baseline_irnt.regenie regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt''')

for chr in range(2,23):
    os.system(f"""awk 'FNR>1 {{print $0}}' regenie_ukb_step2_linear_model2_chr{chr}_fiall_baseline_irnt.regenie >> regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt""")

os.system("gzip regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt")
