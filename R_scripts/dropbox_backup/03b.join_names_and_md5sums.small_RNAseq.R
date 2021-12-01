### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Documentation")

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

# file names path
file_names_path <- file.path(inpath, "file_names.testis.small_RNAseq.csv")

# md5 sums path
md5_sums_path <- file.path(inpath, "file_md5sums.testis.small_RNAseq.csv")

######################################################## READ DATA
# read file names table
file_tb <- readr::read_csv(file_names_path)

# read md5 sums
md5_sums <- readr::read_csv(md5_sums_path, col_names = c("md5_sum", "sample_id"))

######################################################## MAIN CODE
# join, save
file_md5 <- 
  file_tb %>% 
  dplyr::left_join(., md5_sums, by = c("file_name" = "sample_id")) %T>% 
  readr::write_csv(., file.path(outpath, "file_names.md5_sums.testis.small_RNAseq.csv"))







