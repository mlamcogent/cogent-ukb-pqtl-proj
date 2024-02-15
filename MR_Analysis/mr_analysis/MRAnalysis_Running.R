'''
This R script will perform Mendelian randomization analysis 
    using TwoSampleMR and MendelianRandomization

    Also perform Heterogeneity, Directionality and  H pleiotropy tests using TwoSampleMR
'''

suppressWarnings(suppressPackageStartupMessages({
    library(tidyr)
    library(glue)
    library(dplyr)
    library(plyr)
    library(data.table)
    library(magrittr)
    library(TwoSampleMR)
    library(ieugwasr)
    library(mrpipeline)
    library(gwasvcf)
    library(gwasglue)
    library(MRInstruments)
    library(VariantAnnotation)
    library(argparse)
    library(MendelianRandomization)
}))


parser <- ArgumentParser()

parser$add_argument("--exposure_file",
                        help="LD pruned exposure file: eg All_significant_cis_exposure_AfterQC_LDclumping.csv")
parser$add_argument("--output_prefix",
                        help="Output file name; eg: Decode_cis_exposure_PGC3_SCZ")
parser$add_argument("--outcome_file",
                        help="Outcome GWAS file in vcf format; eg PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz")

parser$add_argument("--outcome_fileindex",
                        help="Outcome GWAS file index ; eg PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.index")


args <- parser$parse_args()
exposure_file<-args$exposure_file
prefix=args$output_prefix
outcomevcf=args$outcome_file
outcomevcf_index=args$outcome_fileindex


dbfile="Referencefile/ukb_imp_chr1_22_mac50_info07_b0_7_patched_bfiles_Hg38/ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles"
set_bcftools('Softwares/anaconda3/bin/bcftools')
set_plink('Softwares/Plink/plink')

#method for TwoSampleMR
methodlist_1<-c("mr_wald_ratio","mr_ivw","mr_ivw_mre","mr_ivw_fe")


##Format Exposure
format_exposure <- function(selected_df) {
            exp_df <- format_data(
                selected_df,
                type = "exposure",snps = NULL,header = TRUE,
                phenotype_col = "exposure",snp_col = "ID",
                beta_col = "ES",se_col = "SE",
                eaf_col = "AF",
                effect_allele_col = "ALT",other_allele_col = "REF", # nolint: line_length_linter.
                pval_col = "P",
                samplesize_col = "SS",
                gene_col = "id",min_pval = 1e-800,
                chr_col = "seqnames",pos_col = "start",
                log_pval = FALSE
            )
            return(exp_df) }


#format Outcome data
format_outcome<-function(outcome_file,outcome_fileindex,exposure_df,ldfile){
                        vcf <- query_gwas(outcome_file,
                                        rsid=exposure_df$SNP,
                                        rsidx=outcome_fileindex,
                                        proxies="yes",
                                        bfile=ldfile,
                                        tag_r2=0.8)
                        outcome_df<-gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
                        return(outcome_df)}


###MR Analysis
perform_mr_analysis <- function(dat) {
        mr_res <- data.frame()
        mr_het <- data.frame()
        mr_pleo <- data.frame()
        mr_direction <- data.frame()

    #Two sample MR
        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_1)
            mr_res <- generate_odds_ratios(mr_res)
            mr_res_df <<- rbind.fill(mr_res_df, mr_res)
        }, error = function(e) {
            # Handle the error or perform any desired action
        })

    #heterogeneity
        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogenoty_list_1)
        }, error = function(e) {})

        if (nrow(mr_het)>0) {
            mr_het_df <<- rbind.fill(mr_het_df, mr_het)
        }

    #H pleiotropy
       tryCatch({
            mr_pleo <- mr_pleiotropy_test(dat)
        }, error = function(e) {})

        if (nrow((mr_pleo)>0)) {
            mr_pleo_df <<- rbind.fill(mr_pleo_df, mr_pleo)
        }

    #Directionality
        tryCatch({
            mr_direction <- directionality_test(dat)
        }, error = function(e) {})

        if (nrow(mr_direction)>0) {
            direction_df <<- rbind.fill(direction_df, mr_direction)
        }
    }


##Create empty data frame
mr_res_df <- data.frame() #To save TwoSampleMR MR analysis Results
mr_het_df <- data.frame() #To save TwoSampleMR heterogeneity analysis Results
mr_pleo_df <- data.frame() #To save TwoSampleMR H pleiotropy analysis Results
direction_df <- data.frame() #To save TwoSampleMR Directionality Results
mr_ivw_delta_df<-data.frame() #To save MendelianRandomization ivw MR analysis Results

###Exposure file
exposure_df <- fread(exposure_file)
if ("X" %in% unique(exposure_df$seqnames)) { exposure_df<-exposure_df[exposure_df$seqnames!="X"] }


for ( Exposure in sort(unique(exposure_df$exposure))  ) {
    
    formated_selected_exposure_df <- data.frame()
    formated_outcomedf <- data.frame()
    dat <- data.frame()
    results <- data.frame()
    mr_res <- data.frame()


    #Select the exosure (Protein/Gene)
    tryCatch({exposure_selected_df <- exposure_df[exposure_df$exposure == as.character(Exposure), ]
             }, error = function(e) {})
    
    #Format exposure
    if (nrow(exposure_selected_df) > 0) 
        {formated_selected_exposure_df <- format_exposure(exposure_selected_df)}
    
    if (nrow(formated_selected_exposure_df) > 0) 
        formated_selected_exposure_df$SNP <- toupper(formated_selected_exposure_df$SNP)}

    ## Select the SNPs from outcome GWAS
    tryCatch(
            {if (nrow(formated_selected_exposure_df) > 0) 
                {formated_selected_exposure_df$exposure<-Exposure
                formated_outcomedf <-  format_outcome(outcomevcf, outcomevcf_index,formated_selected_exposure_df, dbfile)
                formated_outcomedf$SNP <- toupper(formated_outcomedf$SNP)}
            }, error = function(e) {}
        )

    if (nrow(formated_outcomedf) > 0) 
        {dat <- harmonise_data(exposure_dat = formated_selected_exposure_df, outcome_dat = formated_outcomedf,action = 1)}

    #Two Sample MR test https://mrcieu.github.io/TwoSampleMR/
    if (nrow(dat) > 0) 
                {
                    perform_mr_analysis(dat)
                }
    
    #MendelianRandomization IVW test  ;https://cran.r-project.org/web/packages/MendelianRandomization/index.html
        if (nrow(dat) > 0) {
            t <- data.frame()
            tryCatch({
            t <- MendelianRandomization::mr_ivw(dat2[[1]], weights = "delta", distribution = "normal")
            }, error = function(e) {})

            if (typeof(t) == "S4") {
            mr_ivw_delta_df <<- rbind.fill(mr_ivw_delta_df, data.frame(
                Model = t@Model,
                Outcome = t@Outcome,
                Penalized = t@Penalized,
                Estimate = t@Estimate,
                CILower = t@CILower,
                Alpha = t@Alpha,
                SNPs = t@SNPs,
                Heter.Stat = paste(unique(t@Heter.Stat),collapse = ","),
                Exposure = t@Exposure,
                Robust = t@Robust,
                Correlation = t@Correlation,
                StdError = t@StdError,
                CIUpper = t@CIUpper,
                Pvalue = t@Pvalue,
                RSE = t@RSE
            ))
            }
        }

write.csv(harmonised_df, glue('{prefix}_Harmonised_Exposure_Outcome.csv'),row.names = FALSE)
write.csv(mr_res_df, glue('{prefix}_TwoSampleMR_Analysis_Multiple_MR_Test.csv'),row.names = FALSE)
write.csv(mr_het_df, glue('{prefix}_TwoSampleMR_Analysis_HeterogenityTest.csv'),row.names = FALSE)
write.csv(mr_pleo_df, glue('{prefix}_TwoSampleMR_Analysis_Hpleiotropy_Test.csv'),row.names = FALSE)
write.csv(direction_df, glue('{prefix}_TwoSampleMR_Analysis_Directionality_Test.csv'),row.names = FALSE)
write.csv(mr_ivw_delta_df, glue('{prefix}_MendelianRandomization_IVW_Delta_Test.csv'),row.names = FALSE)
