### INFO: 
### DATE: 10. 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Rovers_2015_CellRep_GSE64942/Data/Raw/Cleaned")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)
library(doMC)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Rovers_2015_CellRep_GSE64942/Data/Raw/Cleaned"
fasta_path <- file.path(outpath, "adapters.fa")

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# read fasta
adapter <- 
  readr::read_lines(fasta_path) %>% 
  .[!(stringr::str_detect(., ">"))] %>% 
  stringr::str_replace(., "N{4}", "")
  
######################################################## MAIN CODE
# generate all 4-mers, join with adapter
kmers_adapter <- 
  expand.grid(rep(list(c('A', 'G', 'T', 'C')), 4)) %>% 
  do.call(str_c, .) %>% 
  stringr::str_c(., adapter) %>% 
  magrittr::set_names(., str_c("kmer_", 1:length(.))) %>% 
  Biostrings::DNAStringSet(.) %T>% 
  Biostrings::writeXStringSet(., filepath = file.path(outpath, "adapters_kmers.fa"))

