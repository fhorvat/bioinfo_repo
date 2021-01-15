### INFO: 
### DATE: Wed Apr 15 18:41:51 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Documentation")

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
raw_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Raw"

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped"

# genomes path
genome_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genomes"

# genomes
genomes <- list.dirs(genome_path, full.names = F, recursive = F)
  
######################################################## READ DATA

######################################################## MAIN CODE
# list mapped files
mapped_tb <-
  tibble(bam_path = list.files(path = file.path(mapped_path, genomes), pattern = ".*\\.bam$", full.names = T, recursive = F)) %>% 
  dplyr::mutate(genome = str_extract(bam_path, str_c(genomes, collapse = "|"))) %>% 
  dplyr::filter(!(genome %in% c("cow.bosTau9"))) %>% 
  dplyr::mutate(sample_id = bam_path %>% basename(.) %>% str_remove(., "\\.bam$")) %>% 
  dplyr::select(sample_id, bam_path)

# list raw files
raw_tb <- 
  tibble(bam_path = list.files(path = file.path(raw_path, genomes), pattern = ".*\\.gz$", full.names = T, recursive = F)) %>% 
  dplyr::mutate(genome = str_extract(bam_path, str_c(genomes, collapse = "|"))) %>% 
  dplyr::filter(!(genome %in% c("cow.bosTau9"))) %>% 
  dplyr::mutate(sample_id = bam_path %>% basename(.) %>% str_remove(., "_[1,2]{1}\\.txt\\.gz"), 
                pairing = str_extract(bam_path, "(?<=_)[1,2]{1}(?=\\.txt\\.gz)") %>% str_replace_all(., c("1" = "raw_first_path", "2" = "raw_second_path"))) %>% 
  tidyr::pivot_wider(id_cols = "sample_id", names_from = "pairing", values_from = "bam_path") %>% 
  dplyr::left_join(., mapped_tb, by = "sample_id") %>% 
  dplyr::mutate(stage = "GV") %>% 
  dplyr::select(sample_id, stage, raw_first_path, raw_second_path, bam_path)

# save as full table
readr::write_csv(raw_tb, path = file.path(outpath, "maternal_transcriptomes.all_samples.sampleTable.csv"))
