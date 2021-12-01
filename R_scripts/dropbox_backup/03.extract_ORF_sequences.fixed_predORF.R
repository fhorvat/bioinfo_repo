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
consensus_path <- file.path(inpath, "L1_full_length_manual_200707.without_5p_repeat.consensus.fasta")

######################################################## READ DATA
# read consensus sequence
consensus_seq <- Biostrings::readDNAStringSet(consensus_path)

######################################################## MAIN CODE
# # find ORFs
line1_orfs <- predORF(consensus_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = F)

# join to one GRanges object
insertions_orfs_list <- 
  line1_orfs %>% 
  unname(.) %>% 
  do.call(c, .)

# set repeatMasker id
mcols(insertions_orfs_list)$rmsk_id <- as.character(seqnames(insertions_orfs_list))

# split by name and reading frame
insertions_orfs_list <- split(insertions_orfs_list, list(insertions_orfs_list$rmsk_id, insertions_orfs_list$inframe2end))


### reduce ORF coordinates per reading frame and strand in each insertion
# loop
insertions_orfs_reduced <- purrr::map(names(insertions_orfs_list), function(orf_name){
  
  # get one insertion/frame
  insertions_orf <- insertions_orfs_list[[orf_name]]
  
  # get name and reading frame
  rmsk_id <- unique(mcols(insertions_orf)$rmsk_id)
  
  # to GRanges, reduce
  insertions_orf <- 
    insertions_orfs_list[[orf_name]] %>% 
    GenomicRanges::reduce(., ignore.strand = F, min.gapwidth = 0)
  
  # set rmsk_id and reading frames
  mcols(insertions_orf)$rmsk_id <- rmsk_id
  
  # return 
  return(insertions_orf)
  
})

### find two longest ORFs per rmsk_id and save
# get tidy table
insertions_orfs_tb <- 
  insertions_orfs_reduced %>% 
  do.call(c,.) %>% 
  as_tibble(.) %>% 
  arrange(desc(width)) %>% 
  dplyr::slice(1:2) %>% 
  dplyr::summarise(seqnames = unique(seqnames), 
                   start = min(start), 
                   end = max(end), 
                   strand = unique(strand)) 

# get ORFs sequence 
consensus_ORFs <- getSeq(consensus_seq, GRanges(insertions_orfs_tb))
names(consensus_ORFs) <- "L1_full_length_manual_200707.consensus.ORFs"

# save
Biostrings::writeXStringSet(consensus_ORFs, file.path(outpath, "L1_full_length_manual_200707.consensus.ORFs.fasta"))

