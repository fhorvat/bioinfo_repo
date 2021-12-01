### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: Fri May 11 11:36:46 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/documentation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# sample table
sample_table_path <- file.path(outpath, "CNOT6L\ sample\ list\ 11919R_2015-10-29.csv")

# raw files path
raw_path <- "/common/RAW/Svoboda/HamsterUtah/11919R/Fastq"

# mapped files path
mapped_path_mouse <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10"
  
# mapped files path
mapped_path_hamster <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1_vfranke"

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# bam mapped table
bam_table <- 
  tibble(mapped_path = list.files(path = c(mapped_path, mapped_path_hamster), pattern = "*.total.bam$", full.names = T)) %>%
  dplyr::mutate(ID = basename(mapped_path) %>% stringr::str_remove(., ".PE.total.bam"))

# fastq raw table
raw_table <- 
  tibble(raw_path = list.files(path = raw_path, pattern = "*txt.gz$", full.names = T)) %>% 
  dplyr::mutate(pair = ifelse(str_detect(raw_path, "_1.txt.gz$"), "pair_1", "pair_2"), 
                seqname = basename(raw_path) %>% stringr::str_remove(., "_[1,2].txt.gz")) %>% 
  tidyr::spread(data = ., key = pair, value = raw_path) %>% 
  dplyr::mutate(ID = stringr::str_replace(seqname, "_.*", "")) %>% 
  dplyr::left_join(sample_table %>% dplyr::select(ID, stage = `Time Course`, genotype = `Treatment/Control`), by = "ID") %>%
  dplyr::mutate(order = str_extract(ID, "(?<=X)[0-9]+$") %>% as.numeric(.), 
                sample_name = str_c("s_", stage, "_", genotype)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::group_by(sample_name) %>% 
  dplyr::mutate(sample_name2 = str_c(sample_name, "_r", 1:n())) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(ID = sample_name2, stage, genotype, pair_1, pair_2) %>% 
  dplyr::left_join(., bam_table, by = "ID") %T>% 
  readr::write_csv(., file.path(outpath, "path_table.CNOT6L.csv"))

