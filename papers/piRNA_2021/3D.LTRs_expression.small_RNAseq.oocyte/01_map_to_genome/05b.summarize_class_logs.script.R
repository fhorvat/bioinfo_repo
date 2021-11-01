#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## READ DATA
# list files
log_class_list <- list.files(path = outpath, pattern = "s_.*\\.read_stats\\.txt", full.names = T)

######################################################## MAIN CODE
# reads class
log_class_df <-
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>%
  dplyr::bind_rows(.) %>% 
  tidyr::pivot_wider(sample_id, names_from = read_group, values_from = count) %>% 
  dplyr::mutate(mapped_total = rowSums(.[, 2:12])) %>% 
  dplyr::select(sample_id, mapped_total, everything()) %T>% 
  readr::write_delim(x = ., file = "log.read_stats.txt", delim = "\t")
