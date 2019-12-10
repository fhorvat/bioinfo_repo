### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression/Documentation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of experiments
experiment_list <- c("Deng_2014_Science_GSE45719", 
                     "Hamazaki_2015_Development_PRJDB2994", 
                     "Smallwood_2011_NatGenet_PRJEB2547", 
                     "Yamaguchi_2013_CellRes_GSE41908", 
                     "Gan_2013_NatCommun_GSE35005", 
                     "ENCODE_2014_Nature_GSE49417", 
                     "Fugaku",
                     "Veselovska_2015_GenomeBiol_GSE70116")

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq")

# experiment paths
mapped_path <- c(file.path(base_path, experiment_list, "Data/Mapped/STAR_mm10"), 
                 "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_new")
# bam files
bam_path <- list.files(mapped_path, 
                       pattern = ".*\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$|\\.total\\.bam$|\\.bam$", 
                       full.names = T)

# stats and tracks tables
stats_path <- list.files(mapped_path, 
                         pattern = "log\\..*\\.stats_and_tracks\\.csv", 
                         full.names = T)

######################################################## READ DATA
# read stats table
stats_tb <- purrr::map(stats_path, function(path){
  
  # read and tidy stats and tracks
  readr::read_csv(path) %>% 
    dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
    mutate(experiment = str_extract(path, str_c(experiment_list, collapse = "|")))
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
# prepare sample table
sample_table <- 
  tibble(bam_path = bam_path) %>% 
  mutate(sample_name = basename(bam_path),
         sample_id = str_remove(sample_name, "\\.genome\\.Aligned.*\\.bam$|\\.total\\.bam$|\\.bam$"), 
         experiment = str_extract(bam_path, str_c(experiment_list, collapse = "|")), 
         stage = str_remove_all(sample_id, "^s_|_r[0-9]+\\.[S,P]E.*$")) %>% 
  dplyr::left_join(., stats_tb, by = c("sample_id", "experiment")) %>% 
  dplyr::select(sample_id, stage, library_size, sample_name, experiment, bam_path) %T>%
  write_csv(., "developmental_stages.RNAseq.20190426.sampleTable.csv")

# # split  
# sample_tb_list <- split(sample_tb, sample_tb$experiment)
# 
# # write separately 
# purrr::map(names(sample_tb_list), function(experiment){
#   write_csv(x = sample_tb_list[[experiment]], path = str_c("developmental_stages", experiment, "sampleTable.csv", sep = "."))
# })


