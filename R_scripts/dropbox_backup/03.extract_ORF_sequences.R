### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/blast_ORFs_L1")

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
library(Biostrings)
# library(BSgenome.Maur.UCSC.Siomi)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# consensus sequence path
consensus_path <- file.path(inpath, "../blast_consensus_L1/L1_full_length_manual_200707.without_5p_repeat.consensus.fasta")

######################################################## READ DATA
# read consensus sequence
consensus_seq <- Biostrings::readDNAStringSet(consensus_path)

######################################################## MAIN CODE
# # find ORFs
line1_orfs <- predORF(consensus_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)

# get tidy table
line1_orfs_tb<- 
  line1_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  arrange(desc(width)) %>% 
  dplyr::slice(1:2) %>% 
  dplyr::summarise(seqnames = unique(seqnames), 
                   start = min(start), 
                   end = max(end), 
                   strand = unique(strand)) 

# get ORFs sequence 
consensus_ORFs <- getSeq(consensus_seq, GRanges(line1_orfs_tb))
names(consensus_ORFs) <- "L1_full_length_manual_200707.consensus.ORFs"

# save
Biostrings::writeXStringSet(consensus_ORFs, file.path(outpath, "L1_full_length_manual_200707.consensus.ORFs.fasta"))