### INFO: R Script
### DATE: 26.03.2017. 
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################# WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/1cell_motif_discovery")

################################################################# LIBRARIES
library(data.table)
library(stringr)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)

################################################################# PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/1cell_motif_discovery/documentation"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/1cell_motif_discovery/fasta_sequences"
dir.create(outpath, showWarnings = FALSE)

################################################################# SOURCE FILES

################################################################# FUNCTIONS

################################################################# SCRIPT PARAMS

################################################################# TABLES
# read coordinates of 5'UTRs of top 100 transcripts, get GRanges, split by transcriptID
top_5utr <- 
  read_delim(file.path(inpath, "top100_KOvsWT_upreg_ensembl_genes_5UTR_coordinates.txt"), delim = " ", col_names = T) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) 

################################################################################### MAIN CODE
# get sequences, paste together sequences of 5'UTR which have more than one 5'UTR exons, write as FASTA
top_5utr_seq <- 
  split(top_5utr, top_5utr$transcript_id) %>% 
  getSeq(x = Mmusculus, .) %>% 
  lapply(., function(x) as.character(do.call(c, x))) %T>% 
  seqinr::write.fasta(sequences = .,
                      nbchar = 80,
                      names = names(.),
                      as.string = TRUE,
                      file.out = file.path(outpath, "top100_KOvsWT_upreg_ensembl_genes_5UTR_seq.fasta"),
                      open = "w")


