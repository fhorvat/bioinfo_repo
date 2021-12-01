### INFO: 
### DATE: Wed Nov 27 15:10:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/DicerX_embryos/demultiplexed_reads"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path 1
sample_table_raw_path <- file.path(inpath, "Taborska_smalRNA_libs_DicerXXE_191217.csv")

######################################################## READ DATA
# read raw sample data 1
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# clean sample table 1
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::rename_all(str_replace_all, " ", "_") %>% 
  dplyr::select(sample_file_name = Sample_ID, sample_name = `sample_ID_(label)`, genotype, stage = embryo_stage, barcode_seq = Barcode_seq, barcode_id = Barcode_ID) %>% 
  dplyr::filter(!is.na(sample_file_name)) %>% 
  dplyr::mutate(barcode_id = str_c("BC", barcode_id))
                  
# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(barcode_id = raw_path %>% basename(.) %>% str_remove_all(., "_.*")) %>% 
  dplyr::left_join(., sample_table_tidy, by = "barcode_id") %>% 
  dplyr::group_by(genotype, stage) %>% 
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::select(sample_file_name, sample_name, genotype, stage, replicate, barcode_id, barcode_seq, raw_path) %>% 
  dplyr::arrange(genotype, stage, replicate) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, "DicerX", genotype, sample_name, sep = "_") %>% str_c(., "_r", replicate, ".SE"))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, sample_name, barcode_id, barcode_seq, sample_file_name, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c("Taborska_smalRNA_libs_DicerXXE_191217", ".sampleTable.csv"))) %T>% 
  readr:::write_csv(., path = file.path(raw_path, str_c("Taborska_smalRNA_libs_DicerXXE_191217", ".sampleTable.csv")))
  
