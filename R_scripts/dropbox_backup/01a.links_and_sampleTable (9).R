### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/oocytes.Dicer_MT_HET.2021_Apr/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/rat_oocytes.Dicer_MT_HET.RNAseq.2021_Apr"

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
  dplyr::select(tissue, genotype, 
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  as_tibble(.) %>% 
  dplyr::group_by(tissue, genotype) %>% 
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", tissue, genotype, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(genotype = raw_path %>% 
                  basename(.) %>% 
                  str_remove_all(., "\\.fastq\\.gz") %>% 
                  str_remove(., "_S[0-9].*"), 
                genotype = replace(genotype, genotype == "HET", "MT_HET") %>% str_c("Dicer_", .)) %>% 
  dplyr::left_join(., sample_table_tidy, by = "genotype")

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, tissue, genotype, barcode, barcode_no, raw_path) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Dicer_WT", "Dicer_MT_HET"))) %>% 
  dplyr::arrange(tissue, genotype) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



