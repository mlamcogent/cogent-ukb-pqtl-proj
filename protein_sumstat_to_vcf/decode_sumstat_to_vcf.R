
## This R script is to convert the "Decode Pqtl sumstat file to vcf format using gwas2vcf ; https://github.com/MRCIEU/gwas2vcf

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tools))
suppressMessages(library(stringr))
suppressMessages(library(glue))
suppressMessages(library(argparse))
options(scipen=999)


parser <- ArgumentParser(description = "R script for converting decode plasma proteome summary data  to vcf format")

parser$add_argument("--filename", help = "decode plasma proteome summary data ")
args <- parser$parse_args()

##Efffect allele frequency file, 
eaf=fread("assocvariants.annotated_effectallelefrq.txt",sep="\t")
eaf$effectAlleleFreq<-as.numeric(eaf$effectAlleleFreq)

file=args$filename
df<-data.frame()
outname<-substr(file, 1, nchar(file) - 7)
ID<-str_split(outname,"_")[[1]][3]

df=fread(glue("{file}"),sep="\t")
df=df[,c("Chrom","Pos","effectAllele","otherAllele","Beta","SE","minus_log10_pval","N","Name")]
df$minus_log10_pval <- as.numeric(df$minus_log10_pval)

##Since we are mainly focusing on significant markers we selected only the markers that have -log10 P value >=6 # nolint: line_length_linter.
df<- subset(df, minus_log10_pval >= 6)

#if -log10 >200 , it will reduce to 199
df$minus_log10_pval <- ifelse(df$minus_log10_pval < 200, df$minus_log10_pval, 199)

#eaf file merge with original summary stat file
df <- merge(df, eaf, by = "Name")
df=df[,c("Chrom","Pos","Name","effectAllele","otherAllele","Beta","SE","N","minus_log10_pval","effectAlleleFreq")]

#convert -log10 p value to P value
df$minus_log10_pval <- 10^(-df$minus_log10_pval)
df$Chrom <- gsub("chr", "", df$Chrom)

fwrite(df, file = glue("{outname}.tsv"), sep = "\t",row.names = FALSE)
path=getwd()


##convert summary stat file to vcf
r_command <- glue("python ~/software/gwas2vcf/main.py \\
--out {path}/{outname}.vcf \\
--data {path}/{outname}.tsv \\
--ref ~/software/reference_File/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
--dbsnp ~/software/reference_File/GCF_000001405_cleaed.40_Msplited.vcf.gz \\
--json decode_params.txt --id {ID}")

# Execute the R command using 'system'
system(r_command)
system(glue('tabix -f -p vcf {path}/{outname}.vcf.gz'))