### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Ago2_KO/datasets/2019_Oct/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/Ago2_KO_Shubhangini/191021_A00380_0030_AHFWMWDRXX_sdcdLSWxnjsdAKSW"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "Ago2_KO.20191031.sampleTable.raw.csv")

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_ngs_name = Sample_ID, genotype = Description, barcode = index)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "_.*$")) %>% 
  dplyr::left_join(sample_table_tidy, by = c("sample_file_name" = "sample_ngs_name")) %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_c("Ago2_", genotype)) %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, genotype, sample_file_name, "r", sep = "_") %>% str_c(replicate, ".SE"))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, sample_file_name, barcode) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, "Ago2_KO.20191031.sampleTable.csv"))



