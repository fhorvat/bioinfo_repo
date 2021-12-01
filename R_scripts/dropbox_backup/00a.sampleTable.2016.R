### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2016Nov_sequencing")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# inpath
inpath <- getwd()

# sample table
sample_table_path <- list.files(file.path(inpath, "Data/Documentation"), "*sampleTable.csv", full.names = T)

# raw files
raw_path <- list.files(path = c("/common/RAW/Svoboda/2016-11-23-HT3LFBGXY", "/common/RAW/Svoboda/2016-12-12-HVF5KBGXY"), 
                       pattern = "*.txt.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "SE"

# clean raw path
raw_df <- 
  tibble(raw_path) %>% 
  dplyr::mutate(raw_sample_id = raw_path %>% basename(.) %>% str_remove(., "_sequence.txt.gz"))

# clean sample table
sample_table_clean <- 
  sample_table %>% 
  dplyr::select(raw_sample_id = sample, genotype = Type, stage = Sample, barcode = `most common barcode`) %>% 
  dplyr::mutate(genotype = stringr::str_remove(genotype, " B6"), 
                sample_id = str_c("s", genotype, sep = "_"), 
                resequencing = str_detect(raw_sample_id, "HVF5KBGXY")) %>% 
  dplyr::arrange(genotype) %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::mutate(rep = 1:n()) %>%
  dplyr::mutate(sample_id = str_c(sample_id, "_r", rep),
                sample_id = ifelse(resequencing, str_c(sample_id, ".reseq"), sample_id), 
                sample_id = str_c(sample_id, ".", pairing)) %>%
  dplyr::mutate(sample_id = replace(sample_id, sample_id == "s_Lnc21_r3.reseq.SE", "s_Lnc21_r2.reseq.SE"), 
                sample_id = replace(sample_id, sample_id == "s_Lnc21_r4.reseq.SE", "s_Lnc21_r1.reseq.SE"), 
                sample_id = replace(sample_id, sample_id == "s_WT_r3.reseq.SE", "s_WT_r2.reseq.SE")) %>%
  dplyr::mutate(experiment = "2016Nov") %>% 
  dplyr::select(-rep) %>% 
  dplyr::left_join(., raw_df, by = "raw_sample_id") %T>% 
  readr:::write_csv(., path = str_replace(sample_table_path, "\\.csv", ".clean.csv"))

# # only for PE data
# if(pairing == "PE"){
#   
#   sample_table_clean %<>% 
#     rbind(., .) %>% 
#     dplyr::group_by(run) %>% 
#     dplyr::mutate(sequence = 1:n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(run = str_c(run, "_", sequence), 
#                   sample_id = str_c(sample_id, "_", sequence)) %>% 
#     dplyr::select(-sequence)
#   
# }


# save make links script
sample_table_clean %>% 
  dplyr::arrange(sample_id) %>% 
  dplyr::mutate(sample_id = str_c(outpath, "/Data/Raw/Links/", sample_id, ".txt.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", sample_id)) %$%
  make_links %T>% 
  readr::write_lines(., path = str_c(outpath, "/Data/Raw/Links/make_links.sh"))

