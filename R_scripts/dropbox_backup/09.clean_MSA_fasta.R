### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/05_overlap_with_5pUTR_LINE1s")

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

library(Biostrings)
library(DECIPHER)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LINE1s with complete ORFs path
orfs_path <- file.path(inpath, "L1_full_length_manual_200707.consensus.ORFs_blast_hits.csv")

# LINE1s with long 5pUTRs
utrs_path <- file.path(inpath, "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.fasta")

######################################################## READ DATA
# read ORFs table
orfs_tb <- readr::read_csv(orfs_path)

# read 5pUTRs table
utrs_seq <- Biostrings::readDNAStringSet(utrs_path)

######################################################## MAIN CODE
# create sequences without gaps
utrs_without_gaps <-
  utrs_seq %>% 
  DECIPHER::RemoveGaps(., removeGaps = "all", processors = 1)

# save
Biostrings::writeXStringSet(utrs_without_gaps, file.path(outpath, "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.no_gaps.fasta"))



