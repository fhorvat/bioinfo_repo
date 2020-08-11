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
log_list <- list.files(path = outpath, pattern = "s_.*\\.log.out", full.names = T)

######################################################## MAIN CODE
# reads class
log_class_df <-
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>%
  dplyr::bind_rows(.)

# reads logs from bbmap
log_tb <- purrr::map(log_list, function(path){
  
  # read log
  log_tb <- 
    readr::read_delim(path, delim = "\t", skip = 7) %>% 
    dplyr::select(category = `Read 1 data:      `, count = `num reads `, percent = `pct reads`) %>% 
    dplyr::filter(!is.na(count), !str_detect(category, "Rate")) %>% 
    dplyr::mutate(category = str_remove(category, ":.*"),
                  count = str_trim(count) %>% as.numeric(.), 
                  percent = str_remove_all(percent, "%| ") %>% as.numeric(.) %>% divide_by(., 100)) %>% 
    dplyr::slice(1:3) %>% 
    dplyr::mutate(category = str_replace_all(category, c("mapped" = "mapped_total", 
                                                         "unambiguous" = "genome.uniquely_mapped", 
                                                         "ambiguous" = "genome.mapped_to_multiple_loci"))) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.log\\.out")) %>% 
    tidyr::pivot_wider(id_cols = sample_id, names_from = category, values_from = c(count, percent)) %>% 
    dplyr::mutate(raw_input = round(count_mapped_total / percent_mapped_total), 
                  unmapped_total = raw_input - count_mapped_total, 
                  genome.mapped_total = count_genome.uniquely_mapped + count_genome.mapped_to_multiple_loci, 
                  genome.unmapped_total = raw_input - genome.mapped_total) %>% 
    dplyr::select(sample_id, 
                  raw_input, 
                  mapped_total = count_mapped_total, 
                  unmapped_total, 
                  genome.mapped_total, 
                  genome.uniquely_mapped = count_genome.uniquely_mapped, 
                  genome.mapped_to_multiple_loci = count_genome.mapped_to_multiple_loci, 
                  genome.unmapped_total)
  
  # return
  return(log_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# combine with class of reads mapped to genome
if(nrow(log_class_df) > 0){
  
  log_tb %>%
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
  
  log_tb %>%
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
