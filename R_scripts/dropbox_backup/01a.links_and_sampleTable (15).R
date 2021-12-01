### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/golden_hamster.testis.smallRNAseq.8.5dpp.2020/04_joined_by_barcodes"

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
  basename(.) %>% 
  str_remove(., "\\.RNAseq|\\.small_RNAseq")

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name = ID, genotype = Genotype, age = `Age in days`, barcode)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_barcode = raw_path %>% basename(.) %>% str_remove_all(., "^out_|\\.fastq\\.gz")) %>% 
  dplyr::filter(!is.na(sample_barcode), 
                sample_barcode != "unmatched") 

# check if you have all samples
if(!(all(sample_tb$sample_barcode %in% sample_table_tidy$barcode))){
  
  # warn me
  warning("You lost some samples mate! Check what's wrong!")
  
}

# join by sample name
sample_tb %<>% 
  dplyr::left_join(., sample_table_tidy, by = c("sample_barcode" = "barcode")) %>% 
  dplyr::mutate(stage = "testis",
                genotype = str_c("Mov10l_", genotype)) %>% 
  dplyr::group_by(genotype, age) %>% 
  dplyr::mutate(replicate = 1:n(), 
                sample_id = str_c("s", stage, genotype, age, "dpp", sample_name, str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, age, replicate, barcode = sample_barcode, sample_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



