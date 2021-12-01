### INFO: 
### DATE: Mon Jul 23 09:19:25 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Documentation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# sample table path
sample_table_path <- list.files(inpath, "*.samples.txt", full.names = T)

# raw path
raw_path <- "/common/RAW/Svoboda/old_data/T3T_DcrTrans_2011.SOLiD.smallRNASeq"

# set experiment
experiment <- "T3T_DcrTrans_2011" 

######################################################## READ DATA
# read file list
file_list <- readr::read_delim(sample_table_path, delim = "\t")

######################################################## MAIN CODE
# raw table 
raw_table <- 
  tibble(raw_path = list.files(raw_path, "*.gz$", full.names = T)) %>% 
  dplyr::mutate(raw = basename(raw_path) %>% str_remove(., ".gz$"))

# filter sample table
sample_table_filt <- 
  file_list %>% 
  dplyr::mutate(raw = basename(raw)) %>% 
  dplyr::left_join(., raw_table, by = "raw") %>% 
  dplyr::select(-raw) %>% 
  dplyr::mutate(file_type = str_extract(raw_path, "csfasta|qual"), 
                sample_id = str_remove(sample_id, ".csfasta|.qual")) %>% 
  tidyr::spread(key = file_type, value = raw_path) %>% 
  dplyr::mutate(genotype = 
                  str_remove_all(sample_id, "Utra|_r[1-9]$|Utra|T3T_(?=Dcr)") %>% 
                  str_replace(., "T3T_", "WT"), 
                transfection = 
                  str_extract(sample_id, "Utra") %>% 
                  str_replace(., "Utra", "utran") %>% 
                  replace(., is.na(.), "tran"), 
                cell = str_extract(sample_id, "T3T")) %>% 
  dplyr::group_by(genotype, transfection, cell) %>%
  dplyr::mutate(replicate = 1:n(),
                replicate = str_c("r", replicate)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_c("s", cell, genotype, transfection, replicate, sep = "_") %>% str_c(., ".SE")) %>% 
  dplyr::select(sample_id, genotype:cell, csfasta, qual) %T>%
  readr::write_csv(., path = file.path(outpath, str_c(experiment, ".sample_table.csv")))

# write links script
sample_table_filt %>% 
  dplyr::mutate(csfasta = str_c("ln -s ", csfasta, " ", str_c(sample_id, ".csfasta.gz")), 
                qual = str_c("ln -s ", qual, " ", str_c(sample_id, ".qual.gz"))) %>% 
  tidyr::gather(file, links, -c(sample_id:cell)) %$%
  links %>% 
  readr::write_lines(., path = file.path(outpath, "../Raw/Links", str_c(experiment, ".make_links.sh")))

