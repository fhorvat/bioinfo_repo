### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/Dicer_Mili_KO_Eliska/joined_reads"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "Dicer_Mili_KO.sampleTable.raw.csv")

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)
  
######################################################## MAIN CODE
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_ngs_name = `Sample NGS ID`, sample_basespace_name = `Basespace ID`, barcode = `Barcode seq`) %>% 
  dplyr::mutate(sample_basespace_name = str_remove(sample_basespace_name, "^[0-9]+"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "_S[0-9]+_r[1,2]+\\.txt\\.gz")) %>% 
  dplyr::left_join(sample_table_tidy, by = c("sample_file_name" = "sample_basespace_name")) %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_extract(sample_ngs_name, "DBL|MILI|SOM|WT"), 
                replicate = str_extract(sample_ngs_name, str_c("(?<=", genotype, ") *[1-9]o*")) %>% str_squish(.), 
                age = ifelse(str_detect(replicate, "o"), "_old", ""), 
                replicate = str_remove(replicate, "o"), 
                sample_id = str_c("s", stage, genotype, sep = "_") %>% str_c(., age, "_r", replicate, ".PE"), 
                bc = str_extract(sample_ngs_name, "BC[0-9]+(?= )"), 
                pairing = rep(1:2, times = 12))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, "_", pairing, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, age, replicate, bc, barcode, sample_ngs_name, sample_file_name) %>% 
  dplyr::mutate(age = ifelse(age == "_old", "old", "young")) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c("Dicer_Mili_dbl_KO.20190923", ".sampleTable.csv")))
  
  

