### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_MII_Piwil1.RNAseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/golden_hamster.PiwiL1_KO.MII.Siomi.2020"

# links path 
links_path <- file.path(inpath, "../Raw/Cleaned") 

# list samples
sample_list <- list.files(links_path, "\\.fastq\\.gz$", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.) %>% 
  str_remove(., "\\.RNAseq|\\.small_RNAseq")

# create sample table
links <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.clean\\.fastq\\.gz"), 
                read_in_pair = str_extract(sample_file_name, "(?<=R)1$|(?<=R)2$"), 
                stage = "MII", 
                genotype = str_extract(sample_file_name, "hete|homo") %>% str_replace_all(., c("hete" = "Piwil1_HET", "homo" = "Piwil1_KO")), 
                replicate = str_extract(sample_file_name, "(?<=_MII_)[0-9]{1}")) %>% 
  dplyr::mutate(sample_id = str_c("s_", stage, "_", genotype, "_r", replicate, ".PE", "_", read_in_pair)) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz"))

links %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
links %>% 
  dplyr::select(sample_id, stage, genotype, replicate, sample_file_name) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "_[1,2]$"), 
                sample_file_name = str_remove(sample_file_name, "\\.R[1,2]$")) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



