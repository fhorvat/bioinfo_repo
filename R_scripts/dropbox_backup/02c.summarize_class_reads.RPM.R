#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/small_RNAseq_read_distribution/datasets/hamster_oocyte_Mov10l.smallRNAseq")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# list read class files
read_class_path <- list.files(inpath, "\\.read_class\\.widths\\.csv", full.names = T, recursive = T)

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers/6_filter_18to32nt"
library_size_path <- file.path(library_size_path, "library_sizes.txt")

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

# library sizes table
library_sizes <- 
  readr::read_delim(library_size_path, delim = "\t",col_names = c("sample_id", "library_size")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.18to32nt"))

######################################################## MAIN CODE
# set feature classes
class_hier <- 
  tibble(read_group = c("miRNA",
                        "rRNA", "tRNA", 
                        "SINE", "LINE", "LTR", "other_repeat",
                        "protein_coding", "other_ensembl", 
                        "not_annotated"), 
         class = c("miRNA", 
                   "rRNA", "tRNA", 
                   "repeats", "repeats", "repeats", "repeats", 
                   "mRNA", "other_mapped", 
                   "other_mapped")) %>% 
  dplyr::mutate(class = read_group)

# filter table
read_class_tb <- 
  read_class_list %>% 
  dplyr::filter(str_detect(sample_id, "oocytes")) %>% 
  dplyr::filter(!str_detect(sample_id, "19to32nt")) %>% 
  dplyr::filter(read_width >= 18, read_width <= 32) %>% 
  dplyr::left_join(., class_hier, by = "read_group") %>% 
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_sizes, by = "sample_id") %>% 
  dplyr::mutate(rpm = (count / round(library_size / 1E6, 6)),
                class = factor(class, levels = class_hier$class))

# separate for sample
purrr::map(unique(read_class_tb$sample_id), function(sample_name){
  
  # filter, pivot to wide and save
  read_class_sample <- 
    read_class_tb %>% 
    dplyr::filter(sample_id == sample_name) %>% 
    tidyr::pivot_wider(id_cols = c(sample_id, class), names_from = read_width, values_from = rpm) %>% 
    dplyr::arrange(class) %T>% 
    readr::write_csv(., file.path(outpath, str_c(sample_name, "read_class.18to32nt_reads.RPM.20210513.csv", sep = ".")))
  
  # return
  return(sample_name)
  
}) 
