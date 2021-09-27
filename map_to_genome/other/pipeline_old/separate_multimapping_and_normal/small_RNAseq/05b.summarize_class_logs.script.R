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
log_star_list  <-  list.files(path = outpath, pattern = "*.Log.final.out", full.names = T)

######################################################## MAIN CODE
# reads class
log_class_df <-
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>%
  dplyr::bind_rows(.)

# log from STAR
log_final <-
  suppressMessages(purrr::map(log_star_list, read_STAR.Log.final.out, reshape = T)) %>%
  dplyr::bind_rows(.) %>%
  dplyr::mutate_all(.funs = list(~replace(., is.na(.), 0))) %>%
  dplyr::rename(sample_id = log_id) %>%
  magrittr::set_colnames(stringr::str_remove_all(colnames(.), "_reads|number_of_")) %>%
  dplyr::mutate(mapped_total = uniquely_mapped + mapped_to_multiple_loci,
                genome.mapped_total = mapped_total,
                genome.unmapped_total = unmapped) %>%
  dplyr::select(sample_id,
                raw_input = input, mapped_total, unmapped_total = unmapped,
                genome.mapped_total, genome.uniquely_mapped = uniquely_mapped,
                genome.mapped_to_multiple_loci = mapped_to_multiple_loci, genome.unmapped_total)

# combine with class of reads mapped to genome
if(nrow(log_class_df) > 0){
  
  log_class_df %>% 
    tidyr::pivot_wider(sample_id, names_from = read_group, values_from = count) %>%
    dplyr::mutate(mapped_total = rowSums(.[, 2:12])) %>%
    dplyr::select(sample_id, mapped_total, everything()) %>% 
    dplyr::mutate(sample_original = str_remove(sample_id, "(?<=[P,S]E).[0-9]{2}to[0-9]{2}nt$")) %>% 
    dplyr::left_join(., log_final %>% dplyr::select(sample_original = sample_id, raw_input), by = "sample_original") %>%
    dplyr::mutate(unmapped_total = raw_input - mapped_total) %>%
    dplyr::select(sample_id, raw_input,
                  mapped_total, unmapped_total, miRNA.mature.sense:not_annotated) %>%
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., file = "log.read_stats.txt", delim = "\t")
  
}else{
  
  log_final %>%
    dplyr::select(sample_id, raw_input, mapped_total, unmapped_total) %>% 
    dplyr::mutate(miRNA.mature.sense = NA, miRNA.other.sense = NA,
                  protein_coding.sense = NA, rRNA = NA, 
                  SINE = NA, LINE = NA, LTR = NA, other_repeat = NA, 
                  annotated_pseudogene = NA, other = NA, not_annotated = mapped_total) %>% 
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., file = "log.read_stats.txt", delim = "\t")
  
}

