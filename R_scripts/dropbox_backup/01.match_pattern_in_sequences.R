### INFO: 
### DATE: Fri Feb 23 10:14:56 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/splice_sites")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(Biostrings)
library(msa)

######################################################## PATH VARIABLES
outpath <- getwd()

# sequences path
seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/sequences/mammals_plusStrand"

# info table coordinates
info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/Ago2_mammals.csv"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read sequences
sequences <- Biostrings::readDNAStringSet(filepath = list.files(path = seq_path, pattern = "*.fasta", full.names = T))

# read info table
info_df <- readr::read_csv(file = info_path)

######################################################## MAIN CODE
## filter sequences - take only one intron
intron_filt <- "intron3_4"
sequences_filt <- sequences[str_detect(names(sequences), intron_filt)]

# save only intron3_$


### get reverse complement of sequences of Ago2 which are on minus strand
# get names of species which have Ago2 on minus strand, collapse to pattern
minus_strand_animals <- 
  info_df %>%
  dplyr::select(ensembl_dataset, Ago2_strand) %>% 
  dplyr::filter(str_detect(Ago2_strand, "-")) %$%
  ensembl_dataset %>% 
  stringr::str_c(., collapse = "|")

# get index of sequences on minus strand
minus_strand_index <- str_detect(string = names(sequences_filt), pattern = minus_strand_animals)

# set those sequences to their reverse complement
sequences_filt[minus_strand_index] <- Biostrings::reverseComplement(sequences_filt[minus_strand_index])
# Biostrings::writeXStringSet(x = sequences_filt, filepath = "Ago2.intron3_4.mammals.sense.fasta")


# #### align sequences
# # chose animals to align
# chosen_animals <- 
#   info_df %>%
#   dplyr::select(ensembl_dataset, family) %>% 
#   dplyr::filter(family %in% c("Muridae", "Cricetidae")) %$%
#   ensembl_dataset %>% 
#   stringr::str_c(., collapse = "|")
# 
# # filter by animals
# sequences_filt2 <- sequences_filt[str_detect(string = names(sequences_filt), pattern = chosen_animals)]
# sequences_filt2_aligned <- 
#   msa::msa(inputSeqs = sequences_filt2, method = "ClustalOmega") %>% 
#   unmasked(.) %>% 
#   Biostrings::writeXStringSet(x = ., filepath = "Ago2.intron3_4.Muridae.Cricetidae.ClustalO.msa.fasta")


### find splice pattern in sequences
### set splice pattern
# YTNAY (branch sequence 20-50 nucleotides upstream of acceptor site) 
# Y-rich-NCAG/G (acceptor site) (AAG/GAA)
# TGTAAGY (splice donor)

splice_pattern <- "YTNAY"

# find position of splice site
splice_pos <- Biostrings::vmatchPattern(pattern = splice_pattern, subject = sequences_filt, fixed = F, max.mismatch = 0)

# get sequence matching pattern, filter (remove 0-length sequences and sequences containing N)
splice_seq <- 
  sequences_filt[splice_pos] %>% 
  .[!(stringr::str_detect(as.character(.), "N"))] %>% 
  .[nchar(.) > 0]

# filter positions, get data.frame of relative positions
splice_relPos <- 
  splice_pos[names(splice_pos) %in% names(splice_seq)] %>% 
  unlist(.) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(names = str_replace(names, "_intron.*", "")) %>% 
  dplyr::rename(start_rel = start, end_rel = end, ensembl_dataset = names)

# get absolute genomic positions of splice site
splice_absPos <- 
  info_df %>% 
  dplyr::select(scientific_name, ensembl_dataset, coordinates = dplyr::contains(intron_filt), strand = Ago2_strand) %>% 
  tidyr::separate(col = coordinates, into = c("seqnames", "coordinates"), sep = ":") %>% 
  tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-") %>% 
  dplyr::right_join(splice_relPos, by = "ensembl_dataset") %>% 
  dplyr::mutate_at(.vars = vars(matches("start|end")), .funs = funs(as.integer(.))) %>% 
  mutate_cond(strand == "+", start_abs = start + start_rel - 1, end_abs = start_abs + width - 1) %>% 
  mutate_cond(strand == "-", start_abs = end - end_rel + 1 , end_abs = start_abs + width - 1) 
  


