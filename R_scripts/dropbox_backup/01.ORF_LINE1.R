### INFO: 
### DATE: Thu Nov 07 07:19:02 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

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

library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# rmsk path
rmsk_path <- "C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework/rmsk.mm10.20180919.chr11.fa.out.csv"

######################################################## READ DATA
# read data
rmsk_tb <- readr::read_csv(rmsk_path)

######################################################## MAIN CODE
# get LINE1s from RepeatMasker table
line1_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repFamily == "L1")

# get LINE1s as GRanges
# code line with elementNROWS() function will remove insertions with fragments on opposite strands
line1_gr <- 
  line1_tb %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  split(., .$rmsk_id) %>% 
  range(.) %>% 
  .[elementNROWS(.) == 1] %>% 
  unlist(.)
mcols(line1_gr)$rmsk_id <- names(line1_gr)
names(line1_gr) <- NULL

# get ranges of all repeats, remove Low_complexity and Simple_repeat
rmsk_gr <- 
  rmsk_tb %>% 
  dplyr::filter(!(repClass %in% c("Low_complexity", "Simple_repeat"))) %>% 
  GRanges(.)

# find overlaps between LINE1s and other repeats
overlaps <- findOverlaps(line1_gr, rmsk_gr, ignore.strand = T)

# get ID's of LINE1s and overlapping repeats, keep only rmsk_id for ranges which don't overlap any ranges beside their own
overlap_tb <- 
  tibble(rmsk_id = mcols(line1_gr[queryHits(overlaps)])$rmsk_id, 
         insertion_rmsk_id = mcols(rmsk_gr[subjectHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(insertion_rmsk_id = str_c(unique(insertion_rmsk_id), collapse = "|")) %>%
  dplyr::filter(rmsk_id == insertion_rmsk_id)

# keep only LINE1s without any interuptions
line1_without_interuption <- line1_gr[mcols(line1_gr)$rmsk_id %in% overlap_tb$rmsk_id]

# find overlaps between all LINE1 insertions and uninterupted
overlaps <- findOverlaps(line1_without_interuption, line1_gr, ignore.strand = T, type = "within")

# get rmsk ID's to table
overlap_tb <- 
  tibble(rmsk_id = mcols(line1_without_interuption[queryHits(overlaps)])$rmsk_id, 
         insertion_rmsk_id = mcols(line1_gr[subjectHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(insertion_rmsk_id = str_c(unique(insertion_rmsk_id), collapse = "|")) %>%
  dplyr::filter(rmsk_id != insertion_rmsk_id)

# add info about insertion to metadata
mcols(line1_without_interuption)$insert_type <- ifelse(test = mcols(line1_without_interuption)$rmsk_id %in% overlap_tb$rmsk_id, 
                                                       yes = "within", 
                                                       no = "solo")

# get only full length elements
line1_full_length <- line1_without_interuption[width(line1_without_interuption) >= 5500 & width(line1_without_interuption) <= 6500]

# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Mmusculus.UCSC.mm10, names = line1_full_length)
names(line1_seq) <- mcols(line1_full_length)$rmsk_id


# find ORFs
line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)

# get tidy table
line1_orfs_tb<- 
  line1_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(longest_orf_1 >= 3600, 
                longest_orf_2 >= 1100)

# filter LINE1 list
line1_final <- line1_full_length[mcols(line1_full_length)$rmsk_id %in% line1_orfs_tb$rmsk_id]
Biostrings::writeXStringSet(line1_seq[names(line1_seq) %in% mcols(line1_final)$rmsk_id], 
                            filepath = file.path(inpath, "L1_full_length.withORFS.chr11.fasta"))


#### implement ORF finder
# functions
# subseq DNA string by start position
seqReadingFrames <- function(dna_seq, position){
  
  subseq(sequnc, start = position)
  
}

### for each sequence find 2 longest ORFs
seq_top2_orfs <- purrr::map(names(line1_seq), function(dna_seq_names){
  
  # get sequence
  dna_seq <- line1_seq[[dna_seq_names]]
  
  ### get list of sequence and reverse complement of the same sequence 
  seq_rc_list <- list(dna_seq, 
                      reverseComplement(dna_seq))
  
  
  ### iterate over sequence and 3 positions to find all 3 reading frames for each sequence and it's reverse complement
  seq_rf <- purrr::map(seq_rc_list, function(dna_seq){
    
    purrr::map(1:3, function(position){
      
      seqReadingFrames(dna_seq, position)
      
    })
    
  }) %>% 
    unlist(.)
  
  
  ### iterate over all 6 posible reading frames to find the all open reading frame
  seq_orfs_length <- purrr::map(seq_rf, function(dna_seq){
    
    # translate, to character, split, find protein lengths
    dna_seq %>% 
      translate(.) %>% 
      as.character(.) %>% 
      str_split(., pattern = "\\*") %>% 
      unlist(.) %>% 
      sapply(., nchar, USE.NAMES = F)

  })
  
  # for each sequence find 2 longest open reading frames 
  top2_orfs <- 
    seq_orfs_length %>% 
    unlist(.) %>% 
    sort(., decreasing = T) %>% 
    .[1:2]
  
})




