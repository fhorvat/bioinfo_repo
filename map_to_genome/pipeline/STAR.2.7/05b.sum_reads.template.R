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
  
  log_final %>%
    dplyr::left_join(., log_class_df, by = "sample_id") %>%
    dplyr::select(sample_id, 
                  raw_input:genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, 
                  genome.unmapped_total,
                  total, rRNA, repeats, exon, other) %>%
    dplyr::mutate(rDNA_45S.input = NA, rDNA_45S.mapped_total = NA, rDNA_45S.uniquely_mapped = NA,
                  rDNA_45S.mapped_to_multiple_loci = NA, rDNA_45S.unmapped = NA) %>%
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
  
}else{
  
  log_final %>%
    dplyr::mutate(genome.mapped_minus_rDNA = mapped_total,
                  total = NA, rRNA = NA, repeats = NA, exon = NA, other = NA) %>%
    dplyr::select(sample_id, 
                  raw_input:genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, 
                  genome.unmapped_total,
                  total, rRNA, repeats, exon, other) %>%
    dplyr::mutate(rDNA_45S.input = NA, rDNA_45S.mapped_total = NA, rDNA_45S.uniquely_mapped = NA,
                  rDNA_45S.mapped_to_multiple_loci = NA, rDNA_45S.unmapped = NA) %>%
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
  
}

