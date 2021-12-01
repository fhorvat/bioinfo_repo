#!/home/students/fhorvat/R/bin/Rscript
### INFO: read VCF from lnc5 KOs and WT, find exons with mutations in lnc5 KO which are not mutated in WT
### DATE: 22. 5. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling")

######################################################## LIBRARIES
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)
library(tibble)


######################################################## PATH VARIABLES
inpath_GATK <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/results"
inpath_platypus <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/platypus/results"

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
vcf_list_GATK <- list.files(path = inpath_GATK, pattern = "*csv", full.names = T)
vcf_list_platypus <- list.files(path = inpath_platypus, pattern = "*csv", full.names = T)

######################################################## MAIN CODE
mismatches_GATK <- 
  read_csv(file = vcf_list_GATK[1]) %>% 
  dplyr::mutate(fullName = str_c(seqnames, ":", start, "-", end))

mismatches_platypus <- 
  read_csv(file = vcf_list_platypus[1]) %>% 
  dplyr::mutate(fullName = str_c(seqnames, ":", start, "-", end)) %>% 
  dplyr::filter(filter == "PASS")

mismatches_GATK_filtered <-
  mismatches_GATK %>% 
  dplyr::filter(fullName %in% mismatches_platypus$fullName)

mismatches_platypus_filtered <-
  mismatches_platypus %>% 
  dplyr::filter(fullName %in% mismatches_GATK$fullName)

