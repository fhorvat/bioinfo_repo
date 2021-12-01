### INFO: 
### DATE: Wed Nov 27 15:10:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/smallRNASeq_2018/2018-07-05-CCGH0ANXX/Eliska_mESC_MosIR"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "Eliska_mESC_MosIR.sampleTable.raw.csv")

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_file_name = `Sample name`, genotype, transfection = Transfection, cell = `cell type`, 
                barcode_seq = `index sequence read`, NEBnext_index = `NEBnext index no`) %>% 
  dplyr::mutate(genotype = stringr::str_replace(genotype, ", \\?", "_"), 
                transfection = str_replace(transfection, "no transfection", "no_transfection"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = basename(raw_path) %>% 
                  str_remove_all(., "CCGH0ANXX_sM-ES-i1-12_18s003124-1-1_Malik_lane5|_sequence.txt.gz") %>% 
                  str_replace(., "Mos", "_Mos") %>% 
                  str_replace(., "NC", "_NC")) %>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_file_name") %>% 
  dplyr::group_by(genotype, transfection) %>% 
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::select(sample_file_name, genotype, transfection, cell, replicate, barcode_seq, NEBnext_index, raw_path) %>% 
  dplyr::arrange(genotype, transfection, replicate) %>% 
  dplyr::mutate(sample_id = str_c("s", str_remove(sample_file_name, "[1-3]$"), sep = "_") %>% str_c(., "_r", replicate, ".SE"))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, cell, genotype, transfection, replicate, sample_file_name, barcode_seq, NEBnext_index, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c("Eliska_mESC_MosIR", ".sampleTable.csv"))) %T>% 
  readr:::write_csv(., path = file.path(raw_path, str_c("Eliska_mESC_MosIR", ".sampleTable.csv")))

