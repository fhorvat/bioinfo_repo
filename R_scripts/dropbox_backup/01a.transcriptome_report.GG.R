### INFO: 
### DATE: Mon Sep 02 09:57:40 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/transcriptome_assemblies")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
suppressAll <- function(x) suppressMessages(suppressWarnings(x))

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# find Trinity reports
trinity_path <- list.files(inpath, pattern = "s_.*\\.GG\\.Trinity\\.fasta\\.stats", recursive = T)

# find reads-to-transcriptome counts
r2t_counts_path <- list.files(inpath, pattern = "s_.*\\.GG\\..*r2t\\.read_counts\\.txt", recursive = T)

# find transcripts-to-genome counts
t2g_counts_path <- list.files(inpath, pattern = "s_.*\\.GG\\..*t2g\\.read_counts\\.txt", recursive = T)

# find Busco reports
busco_path <- list.files(inpath, pattern = "short_summary_s_.*\\.GG\\.Trinity\\.busco\\.mammalia\\.txt", recursive = T)

######################################################## READ DATA
# read and clean Trinity report
trinity_list <- purrr::map(trinity_path, function(path){
  
  suppressAll(readr::read_delim(path, delim = "\t", col_names = c("category", "count"))) %>%
    dplyr::mutate_all(~(ifelse(is.na(.), "NA", .))) %>%
    dplyr::mutate(count = ifelse(str_detect(category, "Percent GC"), category, count),
                  count = ifelse(str_detect(category, "Total trinity"), str_c(category, count, sep = " "), count)) %>% 
    dplyr::select(-category) %>% 
    dplyr::filter(count != "NA") %>% 
    tidyr::separate(count, into = c("category", "count"), sep = ": ") %>% 
    dplyr::mutate(count = as.numeric(count)) %>% 
    dplyr::slice(1:11)
  
}) %>% 
  set_names(trinity_path %>% str_remove_all(., "/.*|trinity_|\\.GG"))

# read and clean read counts
r2t_counts <- purrr::map(r2t_counts_path, function(path){
  
  tb <- 
    suppressAll(readr::read_delim(path, delim = "\t", col_names = c("count"))) %>% 
    dplyr::mutate(category = c("Reads mapped to transcriptome", "Total reads")) %>% 
    dplyr::select(2:1) %>% 
    tidyr::spread(category, count) %>% 
    dplyr::mutate(`Read mapping ratio` = `Reads mapped to transcriptome` / `Total reads`) %>% 
    tidyr::gather("category", "count")
  
}) %>% 
  set_names(r2t_counts_path %>% str_remove_all(., "/.*|trinity_|\\.GG"))

# read and clean transcripts counts
t2g_counts <- purrr::map(t2g_counts_path, function(path){
  
  tb <- 
    suppressAll(readr::read_delim(path, delim = "\t", col_names = c("count"))) %>% 
    dplyr::mutate(category = c("Transcripts mapped to genome", "Total transcripts")) %>% 
    dplyr::select(2:1) %>% 
    tidyr::spread(category, count) %>% 
    dplyr::mutate(`Transcripts mapping ratio` = `Transcripts mapped to genome` / `Total transcripts`) %>% 
    tidyr::gather("category", "count")
  
}) %>% 
  set_names(t2g_counts_path %>% str_remove_all(., "/.*|trinity_|\\.GG"))

# read and clean BUSCO report
busco_list <- purrr::map(busco_path, function(path){
  
  tb <- 
    suppressAll(readr::read_delim(path, delim = "\t", col_names = c("tmp", "count", "category"))) %>% 
    dplyr::select(category, count) %>% 
    dplyr::mutate(category = ifelse(str_detect(count, "C:[0-9]+\\.[0-9]+\\%\\[.*"), count, category)) %>% 
    dplyr::filter_all(all_vars(!is.na(.))) %>% 
    dplyr::mutate(category = ifelse(str_detect(category, "C:[0-9]+\\.[0-9]+\\%\\[.*"), "Percentage of BUSCOs found", category), 
                  count = ifelse(str_detect(count, "C:[0-9]+\\.[0-9]+\\%\\[.*"), count %>% str_remove_all(., "C:|\\%\\[.*"), count), 
                  count = as.numeric(count))
  
}) %>% 
  set_names(busco_path %>% str_remove_all(., "/.*|trinity_|\\.GG"))

######################################################## MAIN CODE
# bind all tables for one sample
report_list <- purrr::map(names(trinity_list), function(sample_name){
  
  tb <- 
    bind_rows(trinity_list[[sample_name]], 
                  r2t_counts[[sample_name]], 
                  t2g_counts[[sample_name]], 
                  busco_list[[sample_name]]) %>% 
    dplyr::rename_at(vars("count"), ~(sample_name))
  
}) %>% 
  set_names(names(trinity_list))

# join to one table
report_tb <- 
  purrr::reduce(report_list, left_join, by = "category") %T>%
  write_csv(., file.path(outpath, "assemblies_report.csv"))
