### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/datasets")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# set list of dataset names
dataset_name_list <- c("hamster_testis_Mov10l.8.5dpp.small_RNAseq", 
                       "hamster_testis_Mov10l.small_RNAseq", 
                       "hamster_testis_Mov10l.small_RNAseq.reseq", 
                       "hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq")

# get expression values for different datasets
samples_tb_full <- purrr::map(dataset_name_list, function(dataset_name){
  
  # samples path
  base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", dataset_name)
  mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
  samples_list <- list.files(path = mapped_path, pattern = ".*\\.19to32nt\\.bam$", full.names = T)
  
  # create table
  samples_tb <- tibble(sample_id = basename(samples_list) %>% str_remove(., "\\.19to32nt\\.bam"), experiment = dataset_name)
  
  # return
  return(samples_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# save table
readr::write_csv(samples_tb_full, file = file.path(outpath, str_c("Mov10l1_hamster.datasets", "small_RNAseq", "csv", sep = ".")))

