### INFO: finds gaps in assembly (stretches of N longer than 10)
### DATE: Wed Aug 14 21:20:26 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# path to fasta file
fasta_path <- args$fasta_path
fasta_name <- args$fasta_name

######################################################## READ DATA
# read fasta files
scaffold_list <-
  Biostrings::readDNAStringSet(fasta_path) %>%
  as.character(.) %>%
  as.list(.) %>%
  purrr::map(., DNAString)

######################################################## MAIN CODE
# find gaps (N characters)
gaps_tb <- purrr::map(1:length(scaffold_list), function(n){
  
  # extract scaffold
  scaffold <- scaffold_list[[n]]
  
  # find coordinates of gaps
  gap_coords <-
    scaffold %>%
    Biostrings::maskMotif(., "n") %>%
    gaps(.) %>%
    as(., "Views") %>%
    ranges(.)
  
  # ranges
  if(length(gap_coords) > 0){
    
    gap_coords <-
      GenomicRanges::GRanges(seqnames = names(scaffold_list[n]), ranges = gap_coords) %>%
      as.data.frame(.) %>%
      as_tibble(.) %>%
      dplyr::filter(width >= 10) %>%
      dplyr::mutate_all(~as.character(.))
    
  }else{
    tibble()
  }
  
}) %>%
  dplyr::bind_rows(.)

# create bed from data.frame, save
gaps_bed <-
  gaps_tb %>%
  dplyr::mutate(name = width,
                start = as.integer(start) - 1,
                score = 0,
                strand = ".") %>%
  dplyr::select(seqnames, start, end, name, score, strand) %>% 
  dplyr::mutate(seqnames = str_remove(seqnames, " .*")) %T>%
  readr::write_delim(., file.path(outpath, str_c(fasta_name, ".gaps.bed")), delim = "\t", col_names = F)

