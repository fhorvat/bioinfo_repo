### INFO: 
### DATE: Mon May 10 13:58:01 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

######################################################## READ DATA

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = basename(raw_path) %>% str_remove(., "\\.fastq\\.gz"), 
                genotype = str_extract(sample_file_name, "DicerSOM|DicerKO"), 
                sample_number = str_remove(sample_file_name, str_c("-", genotype, ".*")), 
                sample_number_2 = str_remove(sample_file_name, ".*_"),
                mutation = str_remove(sample_file_name, str_c(sample_number, "-")) %>% 
                  str_remove(., str_c(genotype, "-")) %>% 
                  str_remove(., str_c("_", sample_number_2)) %>% 
                  str_remove(., "myc-His-") %>% 
                  str_replace(., "null-mESCs", "null_mESCs"),
                genotype = str_replace(genotype, "Dicer", "Dicer_")) %>% 
  dplyr::group_by(genotype, mutation) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", genotype, mutation, sample_number, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, genotype, mutation, replicate, sample_number, sample_number_2, sample_file_name, raw_path) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



