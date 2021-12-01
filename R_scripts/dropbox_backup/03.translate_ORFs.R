### INFO: 
### DATE: Thu Nov 07 11:23:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/hackaton")

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
library(msa)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# subseq DNA string by start position
subseqStart <- function(dna_seq, position){
  
  subseq(dna_seq, start = position)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# fasta path
fasta_path <- file.path(inpath, "L1_full_length.withORFS.chr11.fasta")

# rmsk path
rmsk_path <- "C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework/rmsk.mm10.20180919.chr11.fa.out.csv"

######################################################## READ DATA
# read fasta, get reverse complement
seq_raw <- Biostrings::readDNAStringSet(fasta_path)
seq_rc <- reverseComplement(seq_raw)

# read repeatMasker
rmsk_tb <- readr::read_csv(rmsk_path)

######################################################## MAIN CODE
### find repNames of LINE1s
rmsk_LINE1 <- 
  rmsk_tb %>% 
  dplyr::select(repName, rmsk_id) %>% 
  unique(.) 

# get vector of names
rmsk_id_names <- 
  rmsk_LINE1$rmsk_id %>% 
  set_names(., rmsk_LINE1$repName)


### get all 6 reading frames for each sequence (3 for original, 3 for reverse complement)
# iterate over sequences
seq_rf <- lapply(names(seq_raw), function(seq_name){
  
  # iterate over 3 possible reading frames
  seq_rf_list <- lapply(1:3, function(frame){
    
    # subseq in 3 positions
    seq_frames <- subseqStart(seq_raw[[seq_name]], position = frame)
    seq_rc_frames <- subseqStart(seq_rc[[seq_name]], position = frame)
    
    # join to one list and return
    return(list(seq_frames, seq_rc_frames))
    
  }) %>% 
    unlist(.)
  
  # return a DNAStringSet list 
  return(DNAStringSet(seq_rf_list))
  
}) %>% 
  set_names(., names(seq_raw))
  

### translate and get longest protein for each sequence
# iterate over sequences
seq_proteins <- lapply(names(seq_rf), function(seq_name){
  
  # translate
  seq_protein <- Biostrings::translate(seq_rf[[seq_name]])
  
  # get the longest protein
  longest_protein <- 
    seq_protein %>% 
    as.character(.) %>% 
    str_split(., "\\*") %>% 
    unlist(.) %>% 
    .[which.max(nchar(.))]
  
  # return
  return(longest_protein)
    
}) %>% 
  unlist(.) %>% 
  AAStringSet(.)

# set names
names(seq_proteins) <- names(seq_rf)
names(seq_proteins) <- names(rmsk_id_names[match(names(seq_proteins), rmsk_id_names)])
  

### multiple sequence alignment and phylogeny tree
# MSA
seq_msa <- msa::msa(seq_proteins)

# convert to seqinr format
seq_msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")

# calculate distance matrix
seq_dist <- dist.alignment(seq_msa_seqinr, "identity")

# get tree and plot
seq_tree <- nj(seq_dist)
plot(seq_tree)
