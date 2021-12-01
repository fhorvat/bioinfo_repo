#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Movl10.RNAseq/Data/Mapped/STAR_mesAur1")

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
# sample list
sample_list <- 
  log_star_list %>% 
  basename(.) %>% 
  str_remove(., "\\.Log\\.final\\.out") %>% 
  stringr::str_remove(., pattern = "\\.genome.*") %>% 
  unique()

# get all posible combinations of alignments (genome/genome.merged)
alignments_all <- 
  expand.grid(sample_list, ".", c("genome", "genome.merged")) %>% 
  do.call(stringr::str_c, .) %>%
  unique(.) %>% 
  tibble(log_id = .)

# reads class
log_class_df <- 
  suppressMessages(purrr::map(log_class_list, readr::read_delim, delim = "\t")) %>% 
  dplyr::bind_rows(.) 

# log from STAR
# log from STAR
log_star <-
  suppressMessages(purrr::map(log_star_list, read_STAR.Log.final.out, reshape = T)) %>%
  dplyr::bind_rows(.) %>%
  dplyr::mutate(log_id = str_replace(log_id, "\\.SE$", ".SE.genome")) %>% 
  dplyr::full_join(., alignments_all, by = "log_id") %>%
  dplyr::mutate_all(.funs = list(~replace(., is.na(.), 0))) %>%
  dplyr::mutate(alignment = str_extract(log_id, pattern = "genome$|genome.merged$"),
         sample_id = str_remove(log_id, pattern = "\\.genome.*")) %>%
  dplyr::select(-log_id) %>%
  dplyr::select(sample_id, alignment, everything()) %>%
  magrittr::set_colnames(stringr::str_remove_all(colnames(.), "_reads|number_of_"))

# mapped reads
log_mapped <-
  log_star %>%
  dplyr::select(-c(unmapped, input)) %>%
  dplyr::mutate(alignment = str_remove(alignment, "\\.merged")) %>%
  dplyr::group_by(sample_id, alignment) %>%
  dplyr::summarise_all(.funs = sum) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(mapped_total = uniquely_mapped + mapped_to_multiple_loci) %>%
  tidyr::gather(category, number_of_reads, -c(sample_id, alignment)) %>%
  tidyr::unite(temp, alignment, category, sep = ".") %>%
  tidyr::spread(temp, number_of_reads) %>%
  dplyr::mutate(mapped_total = genome.mapped_total)

# unmapped reads
log_unmapped <-
  log_star %>%
  dplyr::select(-c(uniquely_mapped, mapped_to_multiple_loci)) %>%
  tidyr::gather(category, number_of_reads, -c(sample_id, alignment)) %>%
  tidyr::unite(temp, alignment, category, sep = ".") %>%
  tidyr::spread(temp, number_of_reads) %>%
  dplyr::mutate(discarded_while_merging = genome.unmapped - genome.merged.input,
                genome.unmapped_total = discarded_while_merging + genome.merged.unmapped)

# combine mapped and unmapped
log_final <- dplyr::left_join(log_mapped, log_unmapped, by = "sample_id")

# combine with class of reads mapped to genome  
if(nrow(log_class_df) > 0){
  
  log_final %>%
    dplyr::left_join(., log_class_df, by = "sample_id") %>%
    dplyr::mutate(unmapped_total = genome.input - mapped_total) %>%
    dplyr::select(sample_id, raw_input = genome.input, mapped_total, unmapped_total,
                  genome.mapped_total, genome.uniquely_mapped, genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, genome.unmapped_total,
                  total, rRNA, repeats, exon, other) %>%
    dplyr::mutate(rDNA_45S.input = NA, rDNA_45S.mapped_total = NA, rDNA_45S.uniquely_mapped = NA, 
                  rDNA_45S.mapped_to_multiple_loci = NA, rDNA_45S.unmapped = NA) %>% 
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
  
}else{
  
  log_final %>%
    dplyr::mutate(unmapped_total = genome.input - mapped_total,
                  genome.mapped_minus_rDNA = mapped_total, 
                  total = NA, rRNA = NA, repeats = NA, exon = NA, other = NA) %>%
    dplyr::select(sample_id, raw_input = genome.input, mapped_total, unmapped_total,
                  genome.mapped_total, genome.uniquely_mapped, genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, genome.unmapped_total,
                  total, rRNA, repeats, exon, other) %>%
    dplyr::mutate(rDNA_45S.input = NA, rDNA_45S.mapped_total = NA, rDNA_45S.uniquely_mapped = NA, 
                  rDNA_45S.mapped_to_multiple_loci = NA, rDNA_45S.unmapped = NA) %>% 
    dplyr::mutate_all(.funs = list(~as.character(.))) %T>%
    readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
  
}
