### INFO: 
### DATE: Fri Jan 04 15:19:25 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other/tardigrade")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# tardigrade genome path
tard_genome_path <- file.path(inpath, "nHd.2.3.abv500.fna")

# blast results path
blast_results_path <- list.files(path = file.path(inpath, "protein_blast"), pattern = ".*\\.blast\\.txt", full.names = T)

######################################################## READ DATA
# read blast results
blast_results <- purrr::map(blast_results_path, function(path){
  
  # read table and name it
  readr::read_delim(file = path, delim = "\t", col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                             "mismatches", "gap_open", 
                                                             "query_start", "query_end", "subject_start", "subject_end", 
                                                             "e_value", "bit_score")) %>% 
    arrange(-alignment_length)
  
}) %>% 
  magrittr::set_names(., basename(blast_results_path) %>% str_remove(., "\\.blast\\.txt"))

# read genome
tard_genome <- Biostrings::readDNAStringSet(filepath = tard_genome_path)

######################################################## MAIN CODE
### GC content
# transform genome to characters
tard_genome_char <- 
  unlist(tard_genome) %>% 
  as.character(.) %>% 
  str_remove_all(., "N")

# calculate GC content
gc_content <- 
  str_count(string = tard_genome_char, pattern = c("G", "C")) %>% 
  sum(.) %>% 
  `/`(., nchar(tard_genome_char)) %>% 
  `*`(., 100)


### dinucleotide composition
# get all dinucleotides
dincleotides <- 
  expand.grid(rep(list(c('A', 'G', 'T', 'C')), 2)) %>% 
  do.call(str_c, .)

# get dinucleotide composition
dinucleotide_composition <- str_count(string = tard_genome_char, pattern = dincleotides)

# make tibble
dinucleotide_composition_tb <- 
  tibble(dinucleotide = dincleotides, count = dinucleotide_composition) %>% 
  mutate(freq = count / sum(dinucleotide_composition) * 100) %>% 
  select(-count) %T>%
  write_csv(., path = file.path(outpath, "nHd.2.3.abv500.dinucleotide_composition.csv"))


### get 10 random 10kb chunks from genome
# filter 10 scaffolds from genome
tard_genome_filt <- tard_genome[width(tard_genome) >= 10000]
tard_genome_filt <- tard_genome_filt[sample(length(tard_genome_filt), 10)]

# sample width
sample_width <- width(tard_genome_filt) - 10000
sample_width <- sapply(sample_width, function(x) sample(1:x, 1))

# subset genome
tard_genome_filt <- Biostrings::subseq(tard_genome_filt, start = sample_width, end = sample_width + 10000 - 1)
names(tard_genome_filt) <- str_c(names(tard_genome_filt), ":", sample_width, "-", sample_width + 10000 - 1)

# write as fasta
Biostrings::writeXStringSet(tard_genome_filt, filepath = file.path(outpath, "nHd.2.3.abv500.10kb.10_chunks.fasta"))
