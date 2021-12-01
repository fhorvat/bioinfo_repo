### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/mouse_testis.Papd7.small_RNAseq.Apr_2021/Data/Documentation")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw files path
raw_path <- "/common/RAW/Svoboda/Papd7_KO_Zuzka/mouse_testis.Papd7.small_RNAseq.Apr_2021"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw\\.csv", full.names = T)

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.)

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name = sample, genotype, age, 
                barcode = `Barcode seq`, barcode_no = `Barcode number`, 
                samples_order) %>% 
  dplyr::mutate(stage = "testis", 
                genotype = genotype %>% str_replace_all(., " ", "_") %>% str_c("Papd7_", .),
                samples_order = as.character(samples_order)) %>% 
  dplyr::group_by(genotype, age) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", stage, genotype, age, sample_name, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::filter(!str_detect(raw_path, "oocytes")) %>% 
  dplyr::mutate(sample_file_name = basename(raw_path) %>% 
                  str_remove(., "\\.fastq\\.gz$"), 
                samples_order = str_extract(sample_file_name, "^[0-9]+"), 
                genotype = str_extract(sample_file_name, "WT|HET|KO(?=-)|KO4") %>% 
                  str_c("Papd7_", .) %>% 
                  str_replace(., "Papd7_KO4", "Papd7_ex4_KO"), 
                age = str_extract(sample_file_name, "p7(?=_)|p14(?=_)") %>% 
                  str_remove(., "^p") %>% 
                  str_c(., "dpp"))

# join by sample name
sample_tb %<>% 
  dplyr::left_join(., sample_table_tidy, by = c("samples_order", "genotype", "age")) 

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, tissue = stage, genotype, age, replicate, barcode, barcode_no, sample_name, sample_file_name, samples_order) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = "."))) %T>% 
  readr::write_csv(., file = file.path(raw_path, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))


