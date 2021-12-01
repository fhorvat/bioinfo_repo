### INFO: compare expression of repeatMasker elements with max. 10 multimap positions (original) and without limiting multimapping (multimap)
### DATE: 25. 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis"
multimap_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis/FPKM_rptmsk_multimap_reads.csv"
original_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis/FPKM_rptmsk_original_reads.csv"
ltr_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_genomic_insertion_freq/LTR_families_data.csv"
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

######################################################## SOURCE FILES
source(file.path(lib_path, "headt.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# read FPKM expression
original_fpkm <- readr::read_csv(file = original_path)
multimap_fpkm <- readr::read_csv(file = multimap_path)
  
# get LTR family data
LTR_data <- readr::read_csv(file = ltr_path)

######################################################## MAIN CODE
# join multimap and original expression, save 
fpkm_compare <- 
  left_join(original_fpkm, multimap_fpkm, by = "full_pos") %>% 
  dplyr::select(order(str_replace(colnames(.), "_multimap", ""))) %>% 
  tidyr::separate(full_pos, into = c("pos", "strand", "repClass", "repName"), sep = "\\|") %>% 
  dplyr::mutate(pos = str_c(pos, "|", strand)) %>% 
  dplyr::select(-strand) %>% 
  dplyr::left_join(., LTR_data[, c("repName", "LTR_subclass")], by = "repName") %>% 
  dplyr::select(pos, repClass, repName, LTR_subclass, everything()) %T>% 
  readr::write_csv(., path = file.path(outpath, "expression_multimap_original_full.csv")) %>% 
  dplyr::filter(!is.na(LTR_subclass)) %T>%
  readr::write_csv(., path = file.path(outpath, "expression_multimap_original_LTRs.csv"))
  
  

