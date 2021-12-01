#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: runs MSA from DECIPHER package
### DATE: Wed Dec 01 16:24:39 2021
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
library(data.table)
library(purrr)

library(DECIPHER)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
threads <- as.numeric(args$threads)
in_file <- args$in_file
out_file <- args$out_file

######################################################## READ DATA
# read sequence
protein_seq <- Biostrings::readAAStringSet(in_file)

######################################################## MAIN CODE
# MSA with DECIPHER package
msa_protein <-
  DECIPHER::AlignSeqs(myXStringSet = protein_seq, iterations = 100, refinements = 100, processors = threads, verbose = T) %>%
  AAMultipleAlignment(.)

# save
Biostrings::writeXStringSet(msa_protein@unmasked, filepath = out_file)
