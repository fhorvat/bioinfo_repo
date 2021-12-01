### INFO: 
### DATE: Mon Jul 08 21:52:06 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/blast_results")

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

# blast results path
blast_results_path <- list.files(path = inpath, pattern = ".*\\.blastn\\.txt", full.names = T)

######################################################## READ DATA
# read blast results
blast_results <- purrr::map(blast_results_path, function(path){
  
  # read table and name it
  readr::read_delim(file = path, delim = "\t", col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                             "mismatches", "gap_open", 
                                                             "query_start", "query_end", "subject_start", "subject_end", 
                                                             "e_value", "bit_score")) %>% 
    arrange(-alignment_length)
  
}) %>% 
  magrittr::set_names(., basename(blast_results_path) %>% str_remove(., "\\.blastn\\.txt"))

######################################################## MAIN CODE
