### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE/human/polyA/Data/Raw/Files")

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

# metadata path
metadata_path <- file.path(inpath, "metadata.tsv")

# files download list path
download_list_path <- file.path(inpath, "file_data.txt")

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE/human/polyA/Data/Documentation"

# links path
links_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE/human/polyA/Data/Raw/Links"

######################################################## READ DATA
# read metadata table
metadata_tb <- readr::read_tsv(metadata_path)

# read files download list
download_list <- readr::read_lines(download_list_path)

######################################################## MAIN CODE
# clean and filter metadata
metadata_tidy <- 
  metadata_tb %>% 
  dplyr::select(accession = `File accession`, tissue = `Biosample term name`) %>%
  dplyr::group_by(tissue) %>% 
  dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(replicate == 1)

# filter download list, save
download_list %>% 
  .[str_detect(., str_c(metadata_tidy$accession, collapse = "|"))] %T>%
  readr::write_lines(., file = file.path(outpath, "file_data.filt.txt"))


### create sample table and links script
# set whether reads are single or paired end
pairing <- "PE"

# create rename executable
rename_file <- 
  
  metadata_tidy %>% 
  
  dplyr::rename(
    
    run = accession
    
  ) %>% 
  
  dplyr::mutate(
    
    tissue = str_replace_all(tissue, " ", "_")
    
  ) %>% 
  
  dplyr::arrange(run) %>% 
  
  dplyr::mutate(name = str_c("s", 
                             tissue,
                             sep = "_")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>%
  dplyr::select(sample_id = name, 
                tissue,
                run) %T>%
  readr:::write_csv(., file = file.path(documentation_path, str_c("ENCODE.human.polyA", ".sampleTable.csv"))) %T>%
  readr:::write_csv(., file = file.path(outpath, str_c("ENCODE.human.polyA", ".sampleTable.csv")))

# save make links script
rename_file %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(links_path, "/", name, ".bam"),
                run = str_c("../Files", "/", run, ".bam")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))









