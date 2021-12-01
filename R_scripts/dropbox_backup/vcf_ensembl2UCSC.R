#!/home/students/fhorvat/R/bin/Rscript
### INFO: read VCF, change seqnames from ENSEMBL format (1, 2, 3, ...) to UCSC format (chr1, chr2, chr3, ...)
### DATE: 22. 5. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/05_BQSR")

######################################################## LIBRARIES
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)
library(tibble)

library(GenomicFeatures)
library(VariantAnnotation)

######################################################## PATH VARIABLES
vcf_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/05_BQSR/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/05_BQSR"

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
# read .vcf file
vcf <- readVcf(file = vcf_path)

######################################################## MAIN CODE
