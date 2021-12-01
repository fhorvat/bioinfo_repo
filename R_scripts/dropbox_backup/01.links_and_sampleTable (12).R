### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/datasets/Ago2_KO/2020_Dec/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/Ago2_KO_Shubhangini/Dec_2020"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "201202_Ago2_dSA_seq.sampleTable.raw.csv")

# barcodes path
barcodes_list <- list.files(raw_path, ".*\\.barcodes\\.txt", full.names = T)

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

# read barcodes
barcodes_tb <- purrr::map(barcodes_list, function(path){
  
  # read table, set sample_id
  readr::read_delim(path, delim = " ", col_names = c("count", "barcode")) %>% 
    dplyr::mutate(sample_file_name = path %>% basename(.) %>% str_remove(., "\\.barcodes\\.txt"), 
                  count = as.numeric(count))
  
}) %>% 
  bind_rows(.) %>% 
  dplyr::group_by(sample_file_name) %>% 
  dplyr::top_n(., 1, count)

######################################################## MAIN CODE
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_ngs_name = name, genotype, barcode = Barcodes, ID, DOB) %>% 
  dplyr::left_join(., barcodes_tb, by = "barcode") %>% 
  dplyr::select(-count)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz")) %>% 
  dplyr::left_join(sample_table_tidy, by = c("sample_file_name")) %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_c("Ago2_", genotype)) %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, genotype, sample_file_name, "r", sep = "_") %>% str_c(replicate, ".SE"))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, sample_file_name, sample_ngs_name, barcode, ID, DOB, raw_path) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, "201202_Ago2_dSA_seq.sampleTable.csv"))



