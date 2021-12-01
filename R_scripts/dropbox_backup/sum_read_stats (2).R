#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarize read stats
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats and mapped logs
read_stats_list <- args$read_stats_list
mapped_logs_list <- args$mapped_logs_list

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
  dplyr::select(sample_id, input, mapped_total, unmapped, 
                uniquely_mapped, mapped_to_multiple_loci, 
                mapped_minus_rDNA,
                total, rRNA, repeats, exon, other) %>%
  dplyr::rename(input.fragments = input, mapped_total.fragments = mapped_total, unamapped.fragments = unmapped,
                uniquely_mapped.fragments = uniquely_mapped, mapped_to_multiple_loci.fragments = mapped_to_multiple_loci, 
                mapped_minus_rDNA.fragments = mapped_minus_rDNA) %>% 
  dplyr::rename(total.reads = total, rRNA.reads = rRNA, repeats.reads = repeats, 
                exon.reads = exon, other.reads = other) %>% 
  dplyr::mutate_all(.funs = list(~ as.character(.))) %T>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")

