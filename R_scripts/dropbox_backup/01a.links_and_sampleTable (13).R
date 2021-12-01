### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.RNAseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/golden_hamster.testis.9dpp.2020"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw\\.csv", full.names = T)

######################################################## READ DATA
# read raw sample data
sample_table_raw <- data.table::fread(sample_table_raw_path)

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.) %>% 
  str_remove(., "\\.RNAseq|\\.small_RNAseq")

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  as_tibble(.) %>% 
  dplyr::select(sample_name = ID, genotype = Genotype, age = `Age in days`, barcode = barcode) %>% 
  dplyr::mutate(sample_name = str_replace(sample_name, "^SO", "So"), 
                stage = "testis",
                genotype = str_c("Mov10l_", genotype)) %>%
  dplyr::group_by(genotype, age) %>% 
  dplyr::mutate(replicate = 1:n(), 
                sample_id = str_c("s", stage, genotype, age, sample_name, str_c("r", replicate), sep = "_") %>% str_c(., ".PE"))

# create sample table
links <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.txt\\.gz"), 
                sample_name = str_extract(sample_file_name, pattern = str_c(sample_table_tidy$sample_name, collapse = "|"))) %>% 
  dplyr::filter(!is.na(sample_name)) %>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_name") %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_table_tidy %>% 
  dplyr::select(sample_id, stage, genotype, age, replicate, barcode, sample_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



