### INFO: 
### DATE: Wed Nov 27 15:10:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2020_Apr.brain_spleen/Data/Documentation")

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

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw files path
raw_path <- "/common/RAW/Svoboda/DicerX_viral_infection.brain"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "200421_smallRNA_brain_spleen_DicerX.csv")

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path) 

######################################################## MAIN CODE
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name, genotype, tissue, Barcode_seq) %>% 
  dplyr::mutate(genotype = str_replace(genotype, "DcrX \\+\\/-", "DcrX_HET"), 
                genotype = str_replace(genotype, "WT", "DcrX_WT")) %>% 
  dplyr::mutate(infection = ifelse(tissue == "brain", "TBEV", "CVB3"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_name = raw_path %>% 
                  basename(.) %>% 
                  str_remove(., "\\.txt\\.gz") %>% 
                  str_remove(., "_.*") %>% 
                  str_replace(., "-", " ")) %>% 
  dplyr::right_join(., sample_table_tidy, by = "sample_name") %>% 
  dplyr::group_by(genotype, tissue, infection) %>% 
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::select(genotype, tissue, infection, sample_name, replicate, barcode_seq = Barcode_seq, raw_path) %>% 
  dplyr::mutate(sample_name = str_replace(sample_name, "brain ", "Br")) %>% 
  dplyr::arrange(genotype, tissue, infection, replicate) %>% 
  dplyr::mutate(sample_id = str_c("s", tissue, genotype, infection, sample_name, sep = "_") %>% str_c(., "_r", replicate, ".SE")) %>% 
  dplyr::select(sample_id, everything())

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c("200421_smallRNA_brain_spleen_DicerX", ".sampleTable.csv")))

