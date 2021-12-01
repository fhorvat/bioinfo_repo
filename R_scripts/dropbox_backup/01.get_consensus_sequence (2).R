### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/blast_consensus_L1")

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

library(msa)
library(ape)
library(seqinr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MSA fasta path
msa_fasta_path <- file.path(inpath, "L1_full_length_manual_200707.fasta")

######################################################## READ DATA
# read MSA fasta
msa_fasta <- Biostrings::readDNAMultipleAlignment(msa_fasta_path, format = "fasta")

######################################################## MAIN CODE
# get the consensus
msa_consenus <-
  msa_fasta %>% 
  msaConsensusSequence(., type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
  toupper(.) %>%
  Biostrings::DNAStringSet(.)

# set names
names(msa_consenus) <- "L1_full_length_manual_200707.consensus"

# find the last repeat on 5'UTR = "CCTGCAGATCGCCTGGGAACTGCACAGAGCTCCCAATCCCGGCGGGAGCATCAGGTAAACAGGCCAGCTTCCCGGCCTGTGTTCCAGAC"
# last full repeat ends with 1751, last partial repeat starts with 1752
vmatchPattern("CCTGCAGATCGCCTGGGAACTGCACAGAGCTCCCAATCCCGGCGGGAGCA", msa_consenus, max.mismatch = 5)

# get consensus with last repeat, save
consensus_with_last_repeat <- 
  msa_consenus %>% 
  subseq(., 1752, width(.)) %T>% 
  Biostrings::writeXStringSet(., file.path(outpath, str_c("L1_full_length_manual_200707", "with_one_5p_repeat", "consensus.fasta", sep = ".")))

# get consensus without any repeats, save
consensus_without_repeats <- 
  msa_consenus %>% 
  subseq(., 1752 + 89, width(.)) %T>% 
  Biostrings::writeXStringSet(., file.path(outpath, str_c("L1_full_length_manual_200707", "without_5p_repeat", "consensus.fasta", sep = ".")))

