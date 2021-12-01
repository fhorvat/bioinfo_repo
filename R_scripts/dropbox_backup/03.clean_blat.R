### INFO: 
### DATE: Mon Jul 08 21:52:06 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/blat_results")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# blat .psl path
blat_psl_path <- list.files(path = inpath, pattern = ".*\\.score\\.psl", full.names = T)

# blat .bed path (with strand info)
blat_bed_path <- list.files(path = inpath, pattern = ".*\\.mm10\\.lnc1_3prime\\.bed", full.names = T)

######################################################## READ DATA
# read blat .psl
blat_psl <- 
  purrr::map(blat_psl_path, function(path){
    
    # read and set column names of .psl, arrange by score
    suppressMessages(readr::read_delim(file = path, delim = "\t", col_names = c("seqnames", "start", "end", "subject_coords", "score", "perc_identity")) %>% 
                       dplyr::mutate(start = start + 1, 
                                     species = basename(path) %>% str_remove(., "\\.mm10\\.lnc1_3prime\\.score\\.psl"), 
                                     query_coords = str_c(species, ":", seqnames, ":", start, "-", end)))
    
  }) %>% 
  dplyr::bind_rows(.)

# read blat .bed
blat_bed <- 
  purrr::map(blat_bed_path, function(path){
    
    # read bed
    rtracklayer::import.bed(path) %>% 
      as_tibble(.) %>% 
      dplyr::select(seqnames, start, end, width, strand) %>% 
      dplyr::mutate(species = basename(path) %>% str_remove(., "\\.mm10\\.lnc1_3prime\\.bed"), 
                    query_coords = str_c(species, ":", seqnames, ":", start, "-", end)) %>% 
      dplyr::select(query_coords, width, strand)
    
  }) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# add strand info, get top score for each species
blat_results <- 
  blat_psl %>% 
  dplyr::left_join(., blat_bed, by = "query_coords") %>% 
  tidyr::separate(subject_coords, c("subject_name", "subject_seqnames", "subject_coords", "hit_coords"), sep = ":") %>%
  tidyr::separate(hit_coords, c("hit_start", "hit_end"), sep = "-") %>%
  dplyr::mutate(mm10_width = as.numeric(hit_end) - as.numeric(hit_start)) %>% 
  dplyr::select(seqnames, start, end, width, strand, score, perc_identity, mm10_width, species) %>% 
  dplyr::group_by(species) %>% 
  top_n(1, score) %>% 
  dplyr::ungroup(.) %>% 
  split(., .$species)

# # save as .bed
# purrr::map(names(blat_results), function(species){
# 
#   # filter table, transform to GRanges, save as .bed
#   blat_results[[species]] %>%
#     GRanges(.) %T>%
#     rtracklayer::export.bed(., file.path(outpath, str_c(species, ".mm10.lnc1_3prime.top_result.bed")))
# 
# })

# save as table
bind_rows(blat_results) %>% 
  tidyr::unite(coords, seqnames, start, end, sep = " ") %T>%
  readr::write_csv(., file.path(outpath, "mm10.lnc1_3prime.top_results.blat.csv"))


