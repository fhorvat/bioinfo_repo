### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/oocytes_2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set animal
animal <- "chinese_hamster"

# outpath
outpath <- file.path(getwd(), tolower(animal), "Data/Raw/Links")

# inpath
inpath <- getwd()

# raw path
raw_path <- "/common/RAW/Svoboda/Pig_Hamster_2018/2018-04-03-HFYYVBCX2"

# raw files
raw_files <- list.files(raw_path, "*txt.gz", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# create rename script
rename_file <- 
  tibble(raw_path = raw_files) %>%
  dplyr::mutate(raw_name = basename(raw_path)) %>% 
  dplyr::filter(str_detect(raw_name, "ChHam")) %>% 
  dplyr::mutate(sample_id = str_remove_all(raw_name, "_[1-4]{1}_sequence.txt.gz")) %>% 
  dplyr::mutate(replicate = group_indices(., sample_id)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(pair = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_extract(sample_id, "ChHam"), 
                sample_id = str_c("s_", sample_id, "_r", replicate, ".PE_", pair, ".txt.gz")) %>% 
  dplyr::select(sample_id, raw_path) %T>% 
  readr:::write_csv(., path = file.path(outpath, "../../Documentation", str_c("oocytes_2018.", tolower(animal), ".sample_table.csv")))

# save make links script
rename_file %>% 
  dplyr::mutate(out = file.path(outpath, sample_id), 
                make_links = str_c("ln -s ", raw_path, " ", out)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))

