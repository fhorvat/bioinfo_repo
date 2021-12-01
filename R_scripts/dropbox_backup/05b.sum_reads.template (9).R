#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Dicer_Mili_KO/Mapped/mm10/perfect_reads")

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
log_original_path  <-  file.path(outpath, "..", "log.read_stats.txt")

######################################################## MAIN CODE
# read log from original mapping
log_original <- readr::read_delim(log_original_path, delim = "\t")

# reads class
log_class_df <-
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "(?<=[P,S]E)\\..*$"))

# join, replace values with perfect reads values, save
log_final <- 
  log_original %>% 
  dplyr::select(-c(rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA)) %>% 
  dplyr::left_join(., log_class_df, by = "sample_id") %>% 
  dplyr::select_at(.vars = vars(colnames(log_original))) %>% 
  dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
