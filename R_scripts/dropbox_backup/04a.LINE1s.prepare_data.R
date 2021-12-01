### INFO: 
### DATE: Mon Oct 28 22:19:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test")

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

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# rmsk path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")
rmks_joined_path <- file.path(genome_dir, "rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

rmks_joined_tb <- readr::read_delim(rmks_joined_path, delim = "\t")
line1_joined <- 
  rmks_joined_tb %>% 
  dplyr::filter(seqnames == "chr11", repFamily == "L1")

######################################################## MAIN CODE
# get only chr11
rmsk_chr11 <- 
  rmsk_tb %>% 
  dplyr::filter(seqnames == "chr11") %T>% 
  readr::write_csv(., file.path(outpath, "rmsk.mm10.20180919.chr11.fa.out.csv"))

