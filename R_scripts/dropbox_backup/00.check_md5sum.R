### INFO: 
### DATE: Tue Sep 22 14:42:39 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/Svoboda/piRNA.Zuzka/datasets/bisulfite_sequencing/2020_Sep.test_run")

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

# get md5 sums path
unix_md5_path <- file.path(inpath, "md5_unix.txt")
windows_md5_path <- file.path(inpath, "md5_windows.txt")

######################################################## READ DATA
# read md5 sums
unix_md5 <- readr::read_delim(unix_md5_path, delim = "\t")
windows_md5 <- readr::read_delim(windows_md5_path, delim = "\t")

######################################################## MAIN CODE
# join and check 
all_md5 <- 
  unix_md5 %>% 
  dplyr::left_join(., windows_md5, by = "sample_id")

# check
all(all_md5$md5_unix == all_md5$md5_windows)

