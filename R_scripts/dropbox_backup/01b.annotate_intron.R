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

######################################################## PATH VARIABLES
outpath <- getwd()

# sequences path
seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/sequences/Ago2.intron3_4.mammals.sense.fasta"

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
sequences <- Biostrings::readDNAStringSet(filepath = seq_path)

# read info table
info_df <- readr::read_csv(file = info_path)

######################################################## MAIN CODE
## filter sequences - take only one intron
intron_filt <- "intron3_4"
animal_filt <- "mmusculus"
sequences_filt <- sequences[str_detect(names(sequences), animal_filt)]

### find splice pattern in sequences
### set splice pattern
# YTNAY (branch sequence 20-50 nucleotides upstream of acceptor site) 
# Y-rich-NCAG/G (acceptor site) (AAG/GAA)
# TGTAAGY (splice donor)
seq_pattern <- "YTNAY"

# find position of splice site, get data.frame of relative positions
pattern_pos <- 
  Biostrings::vmatchPattern(pattern = seq_pattern, subject = sequences_filt, fixed = F, max.mismatch = 0) %>% 
  unlist(.) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(names = str_replace(names, "_intron.*", "")) %>% 
  dplyr::rename(start_rel = start, end_rel = end, ensembl_dataset = names)

# get absolute genomic positions of splice site
pattern_absPos <- 
  info_df %>% 
  dplyr::select(scientific_name, ensembl_dataset, coordinates = dplyr::contains(intron_filt), strand = Ago2_strand) %>% 
  tidyr::separate(col = coordinates, into = c("seqnames", "coordinates"), sep = ":") %>% 
  tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-") %>% 
  dplyr::right_join(pattern_pos, by = "ensembl_dataset") %>% 
  dplyr::mutate_at(.vars = vars(matches("start|end")), .funs = funs(as.integer(.))) %>% 
  mutate_cond(strand == "+", start_abs = start + start_rel - 1, end_abs = start_abs + width - 1) %>% 
  mutate_cond(strand == "-", start_abs = end - end_rel + 1 , end_abs = start_abs + width - 1) %>% 
  dplyr::mutate(pattern = seq_pattern) %>% 
  dplyr::select(seqnames, start = start_abs, end = end_abs, strand, pattern) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)



