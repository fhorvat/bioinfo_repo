### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/circRNA_detection/datasets/Fugaku/STARChip")

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

# counts table path
counts_tb_path <- list.files(inpath, pattern = "circRNA.*\\.countmatrix", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Data/Mapped/STAR_mm10"

# find stats and tracks table
sample_tb_path <- list.files(mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

######################################################## READ DATA
# read counts table
counts_tb <- readr::read_delim(counts_tb_path, delim = "\t")

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# clean sample table
sample_tb_clean <- 
  sample_tb %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.PE"), 
                library_size = (library_size / 10e5))

# normalize count table, calculate mean value
count_tb_norm <- 
  counts_tb %>% 
  tidyr::pivot_longer(-cRNA, names_to = "sample_id", values_to = "count") %>% 
  dplyr::left_join(., sample_tb_clean, by = "sample_id") %>% 
  dplyr::mutate(CPM = round((count / library_size), 3)) %>% 
  tidyr::pivot_wider(id_cols = cRNA, names_from = "sample_id", values_from = "CPM") %>% 
  dplyr::mutate(cRNA = str_replace_all(cRNA, "-|:", " ")) %>% 
  dplyr::arrange(desc(s_GV.WE))

# save
readr::write_csv(count_tb_norm, file.path(outpath, "Fugaku.circRNA.5reads.1ind.countmatrix.CPM.csv"))



