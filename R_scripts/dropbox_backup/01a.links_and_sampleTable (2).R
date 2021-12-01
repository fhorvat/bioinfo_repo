### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_spleen.MosIR.small_RNAseq.2021_Jan/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/MosIR.spleen.small_RNAseq.2021"

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
  basename(.)

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_file_name = Name, tissue = sample, genotype, 
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  as_tibble(.) %>% 
  dplyr::filter(!is.na(barcode_no)) %>% 
  dplyr::mutate(sample_name = str_remove(sample_file_name, "_.*"), 
                tissue = "spleen", 
                genotype = genotype %>% 
                  str_replace(., "mosIR", "MosIR") %>% 
                  str_replace(., "wtDicer", "DicerWT") %>% 
                  str_replace(., "PKR -/-", "delPKR") %>% 
                  str_remove_all(., " *\\+")) %>% 
  dplyr::group_by(genotype) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", tissue, genotype, sample_name, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% 
                  basename(.) %>% 
                  str_remove_all(., "\\.fastq\\.gz") %>% 
                  str_remove(., "_S[0-9]+$") %>% 
                  str_replace_all(., "-", "_")) %>% 
  dplyr::filter(!is.na(sample_file_name)) %>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_file_name") %>% 
  dplyr::filter(sample_file_name != "Undetermined")

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, tissue, genotype, replicate, sample_name, sample_file_name, barcode, barcode_no, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



