### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/ovomucin_KO/datasets/ovomucin_KO.2020_Aug/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse.Ln5.Ovo.2020/joined_lanes/ovomucin"

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
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fq\\.gz", full.names = T)) %>% 
  dplyr::mutate(file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fq\\.gz")) %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_extract(file_name, "WT|KO"), 
                genotype = str_c("Ovomucin_", genotype), 
                replicate = str_extract(file_name, "rep0[1-4]") %>% str_remove(., "rep0")) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, genotype, str_c("r", replicate), sep = "_") %>% str_c(., ".SE")) %>% 
  dplyr::filter(!is.na(sample_id))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, file_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



