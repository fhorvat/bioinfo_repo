### INFO: 
### DATE: Mon May 10 13:58:01 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/haematooncology.Sladja/datasets/mouse_cells.CMO_WT.ATACseq.2021_Mar/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse_cells.CMO_WT.ATACseq.2021_Mar"

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
  dplyr::mutate(sample_file_name = basename(raw_path) %>% str_remove(., "\\.fastq\\.gz")) %>% 
  dplyr::mutate(genotype = str_extract(sample_file_name, "CMO|WT"), 
                cell = str_extract(sample_file_name, "LT|ST"), 
                sample_number = str_extract(sample_file_name, "(?<=)S[0-9]{1}(?=_)"), 
                replicate = "r1",
                sequence = str_extract(sample_file_name, "(?<=_)R[0-9]{1}(?=_)") %>% str_remove(., "R")) %>% 
  dplyr::mutate(sample_id = str_c("s", genotype, cell, replicate, sep = "_") %>% str_c(., ".PE")) %>% 
  dplyr::mutate(name = str_c(sample_id, "_", sequence))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", name, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, genotype, cell, replicate, sample_number) %>% 
  dplyr::arrange(sample_id) %>% 
  dplyr::distinct(.) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



