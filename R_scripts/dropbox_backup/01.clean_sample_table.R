### INFO: 
### DATE: Tue Jul 24 14:48:33 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/T3T_mESC_MosIR.smallRNAseq_2016/Data/Documentation")

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
sample_table_path <- list.files(inpath, "*.sample_table.raw.csv", full.names = T)

# raw path
raw_path <- "/common/RAW/Svoboda/2016-12-16-H33J5BBXX"

# set experiment
experiment <- "T3T_mESC_MosIR.smallRNAseq_2016" 

######################################################## READ DATA
# read file list
file_list <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
# raw table 
raw_table <- 
  tibble(raw_path = list.files(raw_path, "*.gz$", full.names = T)) %>% 
  dplyr::mutate(raw = basename(raw_path))

# filter sample table
sample_table_filt <- 
  file_list %>% 
  dplyr::mutate(raw = basename(raw)) %>% 
  dplyr::left_join(., raw_table, by = "raw") %>% 
  dplyr::select(-raw) %>% 
  dplyr::group_by(genotype, transfection, cell) %>%
  dplyr::mutate(replicate = 1:n(),
                replicate = str_c("r", replicate)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_c("s", cell, genotype, transfection, replicate, sep = "_") %>% str_c(., ".SE")) %>% 
  dplyr::select(sample_id, genotype:cell, raw_path) %T>%
  readr::write_csv(., path = file.path(outpath, str_c(experiment, ".sample_table.csv")))

# write links script
sample_table_filt %>% 
  dplyr::mutate(links = str_c("ln -s ", raw_path, " ", str_c(sample_id, ".txt.gz"))) %$%
  links %>% 
  readr::write_lines(., path = file.path(outpath, "../Raw/Links", str_c(experiment, ".make_links.sh")))

