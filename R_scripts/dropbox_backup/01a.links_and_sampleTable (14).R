### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq/Data/Documentation")

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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw files path
raw_path <- "/common/RAW/Svoboda/golden_hamster.testis.8.5dpp.run_2"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

######################################################## READ DATA

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.) %>% 
  str_remove(., "\\.RNAseq|\\.small_RNAseq")

# create sample table
sample_table_tidy <- 
  tibble(raw_path = list.files(raw_path, pattern = "*\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\_r1.fastq\\.gz")) %>% 
  dplyr::select(sample_file_name) %>% 
  unique(.) %>% 
  dplyr::mutate(sample_name = str_extract(sample_file_name, "So[0-9]{3}-M[0-9]+"), 
                sample_name = ifelse(str_detect(sample_file_name, "half"), str_c(sample_name, "_half"), sample_name), 
                stage = "testis", 
                genotype = str_extract(sample_file_name, "WT|KO") %>% str_c("Mov10l1_", .), 
                age = "8.5dpp") %>% 
  dplyr::group_by(genotype) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", stage, genotype, age, sample_name, str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))
  
# create sample table
links <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz"), 
                sample_name = str_extract(sample_file_name, "So[0-9]{3}-M[0-9]+"), 
                sample_name = ifelse(str_detect(sample_file_name, "half"), str_c(sample_name, "_half"), sample_name)) %>% 
  dplyr::filter(!is.na(sample_name)) %>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_name") %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_table_tidy %>% 
  dplyr::select(sample_id, stage, genotype, age, replicate, sample_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



