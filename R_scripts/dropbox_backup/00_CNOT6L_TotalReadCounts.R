#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 16.03.2017
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/")

################################################################################### LIBRARIES
library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(magrittr)

################################################################################### PATH VARIABLES
inpaths <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/documentation"

################################################################################### 
# CNOT6L experiment table
sample_table <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", col_names = T) %>%
  dplyr::select(ID, stage = `Time Course`, treatment = `Treatment/Control`) %>%
  dplyr::mutate(name = str_c("s", stage, treatment, str_replace(ID, "11919", ""), sep = "_")) 

# library size
library_size_table <- 
  tibble(sample_path = list.files(path = inpaths, pattern = "*.bam$", full.names = T, recursive = T), 
         logs_path = list.files(path = inpaths, pattern = "*Log.final.out$", full.names = T, recursive = T)) %>%
  mutate(ID = str_replace_all(sample_path, "^/.*/|_.*", "")) %>%
  dplyr::rowwise() %>%
  mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size)

# join together, select only 1C knock-out (for now)
sample_table <- 
  dplyr::left_join(sample_table, library_size_table, by = "ID") %>% 
  dplyr::select(ID, name, library_size)

# write as table
write_delim(x = sample_table, path = file.path(outpath, "CNOT6L_library_size.txt"), delim = "\t")