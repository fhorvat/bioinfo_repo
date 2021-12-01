### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/golden_hamster.testis.smallRNAseq.reseq.2020/ZuzkaL-20200806--Ln5-Ovo-GHtestis-demultiplexing"

# links path 
links_path <- file.path(inpath, "../Raw.lanes/Links") 

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
  basename(.) 

# sample table tidy
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name = sample, genotype = `genotype Mov10l`, age, file_name)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fq\\.gz", full.names = T)) %>% 
  dplyr::mutate(file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fq\\.gz"), 
                sample_name = str_extract(file_name, pattern = str_c(sample_table_tidy$sample_name, collapse = "|"))) %>% 
  dplyr::mutate(sample_name = replace(sample_name, is.na(sample_name), "Undetermined"), 
                lane = str_extract(file_name, "L00[1-4]$"), 
                stage = "testis",
                genotype = str_extract(file_name, "WT|KO"), 
                genotype = str_c("Mov10l_", genotype), 
                age = str_extract(file_name, "13dpp|21dpp"), 
                replicate = str_extract(file_name, "rep0[1-4]") %>% str_remove(., "rep0")) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, genotype, age, sample_name, str_c("r", replicate, ".", lane), sep = "_") %>% str_c(., ".reseq.SE"), 
                sample_id = replace(sample_id, is.na(sample_id), str_c(file_name[is.na(sample_id)], ".SE")))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, age, sample_name, replicate, lane, file_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.lanes.csv", sep = ".")))



