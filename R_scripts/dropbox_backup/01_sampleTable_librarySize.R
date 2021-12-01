#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 5. 4. 2017.  
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp")

################################################################################### LIBRARIES
################################################################################### 
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(tibble)

################################################################################### PATH VARIABLES
################################################################################### 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/documentation"

sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/documentation/RNAseq_2016_11_23_sampleTable.csv"
  
################################################################################### SOURCE FILES
################################################################################### 

################################################################################### FUNCTIONS
################################################################################### 

################################################################################### SCRIPT PARAMS
################################################################################### 

################################################################################### TABLES
################################################################################### 
# experiment table
sample_table <- 
  read_csv(sample_table_path) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = str_replace(Treatment.Control, " B6", ""), 
         ID = str_replace(ID, "_16.*", ""), 
         name = str_c(str_replace(ID, ".*_[1-8]_", ""), Treatment.Control, sep = "_"))

# library size
library_size_table <- 
  tibble(sample_path = list.files(path = inpath, pattern = "*.bam$", full.names = T, recursive = T), 
         logs_path = list.files(path = inpath, pattern = "*Log.final.out$", full.names = T, recursive = T)) %>%
  mutate(ID = str_replace_all(sample_path, "^/.*/|_16.*", "")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path)

# join together
sample_table <- 
  dplyr::left_join(sample_table, library_size_table, by = "ID") %>% 
  dplyr::select(ID, name, library_size, sample_path) 

# write as table
write_delim(x = sample_table, path = file.path(outpath, "lncKO_2016_library_size.txt"), delim = "\t")
