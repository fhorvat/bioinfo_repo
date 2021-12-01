### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse_brain.TBEV.small_RNAseq.2021_Feb"

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
  as_tibble(.) %>% 
  dplyr::mutate(tissue = "brain", infection = "TBEV", 
                age = "7wk") %>% 
  dplyr::select(sample_file_name = sample, tissue, genotype = `genotype DicerHAoo`, infection, age, 
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  dplyr::mutate(genotype = str_replace(genotype, "Het", "DicerX_HET")) %>% 
  dplyr::group_by(genotype) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", tissue, genotype, infection, sample_file_name, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% 
                  basename(.) %>% 
                  str_remove_all(., "-.*")) %>% 
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
  dplyr::select(sample_id, tissue, genotype, infection, age, replicate, sample_file_name, barcode, barcode_no, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



