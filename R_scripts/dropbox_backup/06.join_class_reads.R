### INFO: 
### DATE: Mon Apr 27 16:32:28 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/Analysis/2020_paper/small_RNA_reads_annotation")

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

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# datasets path
datasets_path <- file.path(inpath, "datasets")

# experiments
experiments <- c("GarciaLopez_2015_RNA_GSE59254", "Yang_2016_SciAdv_GSE83581", "Tam_2008_Nature_GSE10364")

# class reads path
class_reads_path <- file.path(datasets_path, experiments, "7_class_reads.new")

# library size path
library_size_path <- file.path(datasets_path, experiments, "4_library_size")

# list class reads
class_reads_list <- list.files(class_reads_path, "read_class.*", full.names = T)

# list library sizes
library_size_list <- list.files(library_size_path, "library_sizes.txt", full.names = T)

######################################################## READ DATA
# read class reads
class_reads_tb <- purrr::map(class_reads_list, function(path){
  
  # get experiment
  experiment_name <- str_extract(path, str_c(experiments, collapse = "|"))
  
  # read table
  class_reads <- 
    readr::read_delim(path, delim = "\t") %>% 
    dplyr::mutate(experiment = experiment_name) %>% 
    dplyr::rename(sample_id = sample)
  
  # return
  return(class_reads)
  
}) %>% 
  dplyr::bind_rows(.)

# read library sizes
library_size_tb <- purrr::map(library_size_list, function(path){
  
  # get experiment
  experiment_name <- str_extract(path, str_c(experiments, collapse = "|"))
  
  # read table
  library_size <-
    readr::read_delim(path, delim = "\t", col_names = c("sample_id", "library_size")) %>%
    dplyr::mutate(experiment = experiment_name)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE"))

######################################################## MAIN CODE
# join tables
class_reads_final <- 
  class_reads_tb %>% 
  dplyr::left_join(., library_size_tb, by = c("sample_id", "experiment")) %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "s_|\\.21to23nt"), 
                experiment = str_remove(experiment, "_[0-9]{4}.*"), 
                sample_id = str_c(sample_id, experiment, sep = ".")) %>%
  dplyr::mutate(library_size = (library_size / 10e5), 
                CPM = round((count / library_size), 3)) %>% 
  dplyr::select(read_group, count, CPM, sample_id)

# wide count tables
class_reads_counts <- 
  class_reads_final %>% 
  tidyr::pivot_wider(., id_cols = read_group, names_from = sample_id, values_from = count) %T>%
  readr::write_csv(., file.path(outpath, "read_class.21to23nt.counts.csv"))

# wide CPM tables
class_reads_CPM <- 
  class_reads_final %>% 
  tidyr::pivot_wider(., id_cols = read_group, names_from = sample_id, values_from = CPM) %T>%
  readr::write_csv(., file.path(outpath, "read_class.21to23nt.CPM.csv"))

