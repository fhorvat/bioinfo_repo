### INFO: 
### DATE: Thu Apr 25 16:02:32 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/diffExp/Stein_2015_PLoSGenet_GSE57514")

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# set experiment name
experiment <- "Stein_2015_PLoSGenet_GSE57514"

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10_new")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")
  
# stats and tracks path
stats_path <- list.files(mapped_path, "log\\..*\\.stats_and_tracks\\.csv", full.names = T)

######################################################## READ DATA
# get list of bam files
bam_path <- list.files(path = mapped_path, pattern = "*.total.bam$|*.genome.Aligned.sortedByCoord.out.bam$", full.names = T)

# read stats table
stats_tb <- readr::read_csv(stats_path)

######################################################## MAIN CODE
# get sample table
sample_table <- 
  tibble(bam_path) %>% 
  dplyr::mutate(sample_id = bam_path %>% basename(.) %>% stringr::str_remove(., ".genome.Aligned.sortedByCoord.out.bam|.total.bam"), 
                genotype = str_remove_all(sample_id, "s_|_r[123]{1}|.SE")) %>% 
  dplyr::left_join(., stats_tb %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>% 
  dplyr::select(sample_id, genotype, library_size, bam_path) %T>%
  readr::write_csv(., file.path(documentation_path, str_c(experiment, ".sampleTable.csv")))
  


