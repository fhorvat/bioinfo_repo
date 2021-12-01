### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: Fri May 11 11:36:46 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/SRA/sra/Svoboda/2017_download/RNAseq/ENCODE_2014_Nature_GSE49417")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# sample table
sample_table_path <- file.path(outpath, "metadata.tsv")

# download links
download_links <- file.path(outpath, "links.sub.txt")

######################################################## READ DATA
# read sample table
sample_table <- readr::read_tsv(file = sample_table_path)

# read download links
download_table <- readr::read_delim(file = download_links, delim = "\t", col_names = F)

######################################################## MAIN CODE
# clean download table
download_table %<>% 
  dplyr::rename(path = X1) %>% 
  dplyr::mutate(ID = basename(path) %>% stringr::str_remove(., ".fastq.gz"))

# filter sample_table
sample_table_filter <- 
  sample_table %>% 
  dplyr::filter(`File accession` %in% download_table$ID)
  