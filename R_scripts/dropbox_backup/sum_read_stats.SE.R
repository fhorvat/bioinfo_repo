#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarize read stats
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Srinivasan_2016_NatCommun_GSE75246/Data/Mapped/STAR_mm10")

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of read stats and mapped logs
read_stats_list <- list.files(path = outpath, pattern = "*read_stats.txt", full.names = T)
mapped_logs_list <- list.files(path = outpath, pattern = "*.Log.final.out", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# sample list
sample_list <-
  read_stats_list %>%
  basename(.) %>%
  stringr::str_remove(., pattern = "\\.read_stats\\.txt$")

# reads class
log_class_df <-
  lapply(read_stats_list, readr::read_delim, delim = "\t") %>%
  dplyr::bind_rows(.)

# log from STAR
log_star <-
  lapply(mapped_logs_list, read_STAR.Log.final.out, reshape = T) %>%
  dplyr::bind_rows(.) %>%
  rename(sample_id = log_id) %>%
  magrittr::set_colnames(stringr::str_remove_all(colnames(.), "_reads|number_of_")) %>% 
  mutate(mapped_total = uniquely_mapped + mapped_to_multiple_loci)

# combine mapped and unmapped, combine with class of reads mapped to genome
log_final <-
  log_star %>% 
  dplyr::left_join(., log_class_df, by = "sample_id") %>%
  dplyr::mutate(genome.mapped_total = mapped_total, 
                genome.unmapped_total = unmapped, 
                raw_input_rDNA = NA, total_mapped_rDNA = NA, uniquely_mapped_rDNA = NA, 
                multi_mapped_rDNA = NA, unmapped_rDNA = NA) %>% 
  dplyr::select(sample_id, raw_input = input, mapped_total, unmapped_total = unmapped, 
                genome.mapped_total, genome.uniquely_mapped = uniquely_mapped, genome.mapped_to_multiple_loci = mapped_to_multiple_loci, 
                genome.mapped_minus_rDNA, genome.unmapped_total, 
                total, rRNA, repeats, exon, other, 
                raw_input_rDNA, total_mapped_rDNA, uniquely_mapped_rDNA, multi_mapped_rDNA, unmapped_rDNA) %>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")

