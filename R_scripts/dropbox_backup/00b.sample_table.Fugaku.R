### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Documentation")

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

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq")

# experiment paths
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10_new")

# bam files
bam_path <- list.files(mapped_path, 
                       pattern = ".*\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$|\\.total\\.bam$", 
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
    mutate(experiment = "Fugaku")
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
# prepare sample table
sample_tb <- 
  tibble(bam_path = bam_path) %>% 
  mutate(sample_name = basename(bam_path),
         sample_id = str_remove(sample_name, "\\.genome\\.Aligned.*\\.bam$|\\.total\\.bam$"), 
         experiment = "Fugaku", 
         stage = str_remove_all(sample_id, "^s_|_r[0-9]+\\.[S,P]E.*$|.[S,P]E.*$")) %>% 
  dplyr::left_join(., stats_tb, by = c("sample_id", "experiment")) %>% 
  dplyr::select(sample_id, stage, library_size, sample_name, experiment, bam_path) %T>%
  write_csv(., "Fugaku.20190426.sampleTable.csv")

