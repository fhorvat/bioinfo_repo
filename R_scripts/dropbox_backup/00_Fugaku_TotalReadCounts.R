#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 16.03.2017
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/")

################################################################################### LIBRARIES
library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(magrittr)

################################################################################### PATH VARIABLES
inpaths <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/output/documentation"

################################################################################### 
# library size
sample_table <- 
  tibble(sample_path = list.files(path = inpaths, pattern = "*.bam$", full.names = T, recursive = T), 
         logs_path = list.files(path = inpaths, pattern = "*Log.final.out$", full.names = T, recursive = T)) %>%
  mutate(ID = str_replace_all(sample_path, "^/.*/|.bam", "")) %>%
  rowwise() %>%
  mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path)

# write as table
write_delim(x = sample_table, path = file.path(outpath, "Fugaku_library_size.txt"), delim = "\t")