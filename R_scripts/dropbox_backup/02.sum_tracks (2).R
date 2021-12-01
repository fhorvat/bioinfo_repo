### INFO: 
### DATE: Fri Feb 22 12:50:53 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/datasets")

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

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## READ DATA
# list track files
tracks_path <- 
  list.files(inpath, "log\\..*tracks\\.csv", recursive = T, full.names = T) %>% 
  .[str_detect(., ".*mm10_masked.perfect")]

######################################################## MAIN CODE
# join all, save
tracks_all <- 
  map(tracks_path, read_csv) %>% 
  bind_rows(.) %>% 
  write_csv(., file.path(outpath, "LINE1_full_length.20190222.ZJM.RNAseq_tracks.masked_mm10.csv"))
