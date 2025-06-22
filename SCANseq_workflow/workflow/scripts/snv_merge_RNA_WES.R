#!/usr/bin/env Rscript

args = commandArgs(T)
dir_work = args[1]
input_rna_snv = args[2]
input_wes_snv = args[3]
sample_name = args[4]
wes_sample = args[5]

suppressMessages({
    options(stringsAsFactors = FALSE)
    library(ggplot2)
    library(ggsci)
    library(RColorBrewer)
    library(dplyr)
})

setwd(paste0(dir_work))
###  Loading RNA_SNV
rna_snv <- read.table(gzfile(input_rna_snv), header=T, sep="\t")
rna_snv$pos <- rna_snv$pos+1  # the input_rna_snv is 0-based. vcf format is 1-based.
# add CellName info
rna_snv$CellName <- sample_name
# add WES info
rna_snv$WES_Sample <- wes_sample

###  Loading WES vcf file
# read two times the vcf file, first for the columns names, second for the data
vcf_header <- readLines(input_wes_snv)
vcf_header <- vcf_header[-(grep("#CHROM",vcf_header)+1):-(length(vcf_header))]
vcf_header <- unlist(strsplit(vcf_header[length(vcf_header)],"\t"))
vcf_header[1] <- "CHROM"
wes_snv <- read.table(input_wes_snv, stringsAsFactors = FALSE)
colnames(wes_snv) <- vcf_header
wes_snv <- wes_snv[,c(vcf_header[1:10],wes_sample)]
wes_snv$WES_Sample <- wes_sample
colnames(wes_snv) <- c(vcf_header[1:9],"Normal_Sample","Tumor_Sample","WES_Sample")

###  merge RNA_SNV and WES_SNV
rna_snv$chr_pos_sp <- paste0(rna_snv$chrom,"_", rna_snv$pos,"_", rna_snv$WES_Sample)
wes_snv$chr_pos_sp <- paste0(wes_snv$CHROM,"_", wes_snv$POS,"_", wes_snv$WES_Sample)
df_merge <- merge(
    rna_snv[,c("chrom","pos","ref_base","alt_base","ref_count","alt_count","other_count",
        "total_count","ref_ratio","alt_ratio","CellName","chr_pos_sp")],
    wes_snv, by="chr_pos_sp", all=T)
df_merge$refGene <- gsub("Gene.refGene=","",sapply(strsplit(df_merge$INFO, split=";"), function(x){x[17]}))
df_merge$Normal_GT <- sapply(strsplit(df_merge$Normal_Sample, split=":"), function(x){x[1]})
df_merge$Tumor_GT <- sapply(strsplit(df_merge$Tumor_Sample, split=":"), function(x){x[1]})
df_merge$Normal_AD <- sapply(strsplit(df_merge$Normal_Sample, split=":"), function(x){x[2]})
df_merge$Tumor_AD <- sapply(strsplit(df_merge$Tumor_Sample, split=":"), function(x){x[2]})
df_merge$Normal_AF <- sapply(strsplit(df_merge$Normal_Sample, split=":"), function(x){x[3]})
df_merge$Tumor_AF <- sapply(strsplit(df_merge$Tumor_Sample, split=":"), function(x){x[3]})
df_merge$Normal_DP <- sapply(strsplit(df_merge$Normal_Sample, split=":"), function(x){x[4]})
df_merge$Tumor_DP <- sapply(strsplit(df_merge$Tumor_Sample, split=":"), function(x){x[4]})

write.table(df_merge, gzfile(sprintf("%s__%s_snv_merge_RNA_WES_1-based.tsv.gz",sample_name,wes_sample)), row.names=F, quote=F, sep="\t")
