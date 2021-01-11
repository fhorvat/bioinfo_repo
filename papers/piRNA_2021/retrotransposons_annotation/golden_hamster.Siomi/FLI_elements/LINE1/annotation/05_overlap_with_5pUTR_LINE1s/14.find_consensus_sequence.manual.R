### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(msa)
library(DECIPHER)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# fasta path
fasta_path <- file.path(inpath, "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.dashed.fasta")

######################################################## READ DATA
# read fasta
insertions_seq <- Biostrings::readDNAStringSet(fasta_path)

######################################################## MAIN CODE
# pad with "-" at the end to same length
insertions_seq_msa <- 
  insertions_seq %>% 
  as.character(.) %>% 
  str_pad(., width = max(width(insertions_seq)), side = "right", pad = "-") %>% 
  Biostrings::DNAMultipleAlignment(.)

# get the consensus
seq_consensus <-
  msaConsensusSequence(insertions_seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
  toupper(.) %>%
  Biostrings::DNAStringSet(.)

# names
names(seq_consensus) <- "L1_full_length_manual_complete_5pUTRs.consensus"

# save the consensus
Biostrings::writeXStringSet(seq_consensus, file = file.path(outpath, str_c("L1_full_length_manual_200707.manual_complete_5pUTRs.200708", "consensus", "fasta", sep = ".")))


