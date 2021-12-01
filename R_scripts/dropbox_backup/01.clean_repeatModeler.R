
### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/RepeatModeler/RepeatMasker")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# repeatMasker path
rmsk_path <- list.files(path = inpath, pattern = "rmsk.*raw.fa.out.gz")

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
