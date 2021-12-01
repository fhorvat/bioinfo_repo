#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "read_STAR.Log.final.out.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
outpath <- getwd()

# stats from full mapping
full_stats_path <- list.files(file.path(outpath, ".."), pattern = ".*stats_and_tracks.csv", full.names = T)
  
######################################################## READ DATA
# list files
log_class_list <- list.files(path = outpath, pattern = "*read_stats.txt", full.names = T)

# read stats from full mapping
full_stats <- read_csv(full_stats_path)

######################################################## MAIN CODE
# sample list
sample_list <-
  log_class_list %>%
  basename(.) %>%
  stringr::str_remove(., pattern = "\\.perfect\\.read_stats\\.txt") %>%
  unique(.)

# reads class
log_class_df <-
  lapply(log_class_list, readr::read_delim, delim = "\t") %>%
  dplyr::bind_rows(.) %>% 
  mutate(sample_id = str_remove(sample_id, "\\.perfect"))

# clean full stats
full_stats_clean <- 
  full_stats %>% 
  select(sample_id, raw_input) 
  
# combine mapped and unmapped, combine with class of reads mapped to genome
log_final <-
  dplyr::left_join(log_class_df, full_stats_clean, by = "sample_id") %>%
  dplyr::mutate(mapped_total = total, 
                unmapped_total = raw_input - mapped_total,
                genome.mapped_total = NA, genome.uniquely_mapped = NA, genome.mapped_to_multiple_loci = NA, genome.unmapped_total = NA, 
                raw_input_rDNA = NA, total_mapped_rDNA = NA, uniquely_mapped_rDNA = NA, multi_mapped_rDNA = NA, unmapped_rDNA = NA) %>%
  dplyr::mutate(genome.mapped_minus_rDNA = ifelse(str_detect(sample_id, "\\.PE"), round(genome.mapped_minus_rDNA / 2), genome.mapped_minus_rDNA), 
                mapped_total = ifelse(str_detect(sample_id, "\\.PE"), round(mapped_total / 2), mapped_total), 
                unmapped_total = ifelse(str_detect(sample_id, "\\.PE"), round(unmapped_total / 2), unmapped_total)) %>% 
  dplyr::select(sample_id, raw_input, mapped_total:genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, genome.unmapped_total,
                total, rRNA:other,
                raw_input_rDNA:unmapped_rDNA) %>%
  dplyr::mutate_all(.funs = funs(as.character(.))) %>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")

