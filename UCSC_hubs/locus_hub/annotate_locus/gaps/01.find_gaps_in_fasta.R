### INFO: 
### DATE: Wed Aug 14 21:20:26 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/ovomucin_KO/Analysis/2021_paper/ovomucin_locus/AcoCah/annotate_locus/gaps")

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

# fasta path
fasta_path <- list.files(path = file.path(inpath, ".."), pattern = ".*.fa$", full.names = T)

######################################################## READ DATA
# read fasta files
scaffold_list <- 
  purrr::map(fasta_path, function(path){
    
    # read scaffold sequence, transform to list of DNAStrings
    scaffold_DNAString <- 
      Biostrings::readDNAStringSet(path) %>% 
      as.character(.) %>% 
      as.list(.) %>% 
      purrr::map(., DNAString)
    
    # add species name
    names(scaffold_DNAString) <- str_c((path %>% basename(.) %>% str_remove(., "\\..*")), 
                                       names(scaffold_DNAString), sep = ".")
    
    # return 
    return(scaffold_DNAString)
    
  }) %>% 
  unlist(.)

######################################################## MAIN CODE
# find gaps (N characters)
gaps_tb <- purrr::map(1:length(scaffold_list), function(n){
  
  # extract scaffold
  scaffold <- scaffold_list[[n]]
  
  # find coordinates of gaps
  gap_coords <- 
    scaffold %>% 
    Biostrings::maskMotif(., "N") %>% 
    gaps(.) %>% 
    as(., "Views") %>% 
    ranges(.) %>% 
    GenomicRanges::GRanges(seqnames = names(scaffold_list[n]) %>% str_extract(., "(?<=\\.).*"), 
                           ranges = .) %>% 
    as.data.frame(.) %>% 
    dplyr::filter(width >= 10) %>%
    dplyr::mutate_all(~as.character(.)) %>% 
    dplyr::mutate(species = names(scaffold_list[n]) %>% str_remove(., "\\..*"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  split(., .$species)

# save as .bed
purrr::map(names(gaps_tb), function(species){
  
  # create GRanges from data.frame, save as .bed
  gaps_gr <- 
    gaps_tb[[species]] %>% 
    dplyr::select(-species) %>% 
    dplyr::mutate(name = width, 
                  start = as.integer(start) - 1, 
                  score = 0, 
                  strand = ".") %>% 
    dplyr::select(seqnames, start, end, name, score, strand) %T>%
    readr::write_delim(., file.path(outpath, str_c(species, ".gaps.bed")), delim = "\t", col_names = F)
  
})
