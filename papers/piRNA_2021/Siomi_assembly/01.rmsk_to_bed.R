g### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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

library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly"

# rmsk path
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.clean\\.out\\.fa\\.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")
  
######################################################## MAIN CODE
# table to GRanges
rmsk_gr <- 
  rmsk_tb %>% 
  GRanges(.)

# export
rtracklayer::export.bed(rmsk_gr, rmsk_path %>% str_replace(., "\\.gz", ".bed"))



