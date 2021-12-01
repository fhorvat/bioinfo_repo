#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarizes counts in developmental profile over LTRs and LINE1s from repeatMasker
### DATE: Wed Mar 06 11:05:00 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/read_counts")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)
library(purrr)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeatMasker coordinates path
rmsk_documents_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages"

# repeatMasker coordinates as GRangesList .RDS path 
rmsk_gr_path <- file.path(rmsk_documents_path, "rmsk.L1_and_LTRs.filtered.20190306.removed_genes.GRangesList.RDS")

# repeatMasker table path
rmsk_tb_path <- file.path(rmsk_documents_path, "rmsk.L1_and_LTRs.filtered.20190306.csv")

######################################################## READ DATA
# read repeatMasker GRangesList 
rmsk_gr <- readRDS(rmsk_gr_path)

# read repeatMasker table
rmsk_tb <- read_csv(rmsk_tb_path)

######################################################## MAIN CODE
### get specific LTR and LINE families
# convert to genomicRanges
rmsk_filtered_tb <- 
  rmsk_gr %>% 
  unname(.) %>% 
  unlist(.)

# add names as new column
mcols(rmsk_filtered_tb)$rmsk_id <- names(rmsk_filtered_tb)

# convert to table, join with details about each repeat
rmsk_filtered_tb %<>%
  as.tibble(.) %>% 
  left_join(., rmsk_tb %>% dplyr::select(rmsk_id, repName:repFamily), by = "rmsk_id")


### split for different elements 
# LINE 1
line1_all <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(repFamily == "L1") %>% 
  dplyr::mutate(class = "LINE1_all") %>% 
  GRanges(.) 

# Line1 subfamilies = L1Md_T L1Md_A L1Md_Gf L1Md_F2 L1_Mus2 L1_Mus1 L1Md_F
line1_subfamilies <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(repName %in% c("L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F")) %>% 
  dplyr::mutate(class = repName) %>% 
  GRanges(.) 

# IAP all
iap_all <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(str_detect(repName, "IAP")) %>% 
  dplyr::mutate(class = "IAP_all") %>% 
  GRanges(.) 

# IAPEz-Int with LTRs IAPLTR1_Mm
iap_int_ltr <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(repName %in% c("IAPEy-int", "IAPLTR1_Mm")) %>% 
  dplyr::mutate(class = "IAPEy_int_IAPLTR1_Mm") %>% 
  GRanges(.)

# MuERV-int with MT2 LTRs = MERVL-int + MT2_Mm
mervl_mt2mm <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(repName %in% c("MERVL-int", "MT2_Mm")) %>% 
  dplyr::mutate(class = "MERVL_int_MT2_Mm") %>% 
  GRanges(.)

# MTA all LTRs and ints
mta_int_ltr <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(str_detect(repName, "MTA")) %>% 
  dplyr::mutate(class = "MTA_int_LTR") %>% 
  GRanges(.)

# ORR1A all ints and LTR subgroups
orr1a_int_ltr <- 
  rmsk_filtered_tb %>% 
  dplyr::filter(str_detect(repName, "ORR1A")) %>% 
  dplyr::mutate(class = "ORR1A_int_LTR") %>% 
  GRanges(.)

### bind all and save
# bind rows to one table
retro_all_gr <- 
  c(line1_all, line1_subfamilies, 
    iap_all, iap_int_ltr, 
    mervl_mt2mm, mta_int_ltr, 
    orr1a_int_ltr)

