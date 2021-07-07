### INFO: cleans repeatMasker fa.out table
### DATE: Wed Jul 07 16:03:53 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# repeatMasker path
rmsk_path <- list.files(path = inpath, pattern = "rmsk.*\\.raw\\.fa\\.out\\.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_df <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

######################################################## MAIN CODE
### clean and save repeatMasker
rmsk_df %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
  readr::write_delim(str_replace(rmsk_path, "raw", "clean"), delim = "\t")
