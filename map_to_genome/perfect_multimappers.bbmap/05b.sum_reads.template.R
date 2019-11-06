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
# sample list
sample_list <-
  log_class_list %>%
  basename(.) %>%
  str_remove(., "\\.read_stats\\.txt") %>%
  stringr::str_remove(., pattern = "\\.genome.*") %>%
  unique(.)

# reads class
log_class_df <-
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(raw_input = NA, mapped_total = NA, unmapped_total = NA, 
                genome.mapped_total = NA, genome.uniquely_mapped = NA, genome.mapped_to_multiple_loci = NA, genome.unmapped_total = NA, 
                rDNA_45S.input = NA, rDNA_45S.mapped_total = NA, rDNA_45S.uniquely_mapped = NA,
                rDNA_45S.mapped_to_multiple_loci = NA, rDNA_45S.unmapped = NA) %>% 
  dplyr::select(sample_id, raw_input:genome.mapped_to_multiple_loci, 
                genome.mapped_minus_rDNA, genome.unmapped_total, 
                total, rRNA, repeats, exon, other, 
                rDNA_45S.input:rDNA_45S.unmapped) %>% 
  dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
