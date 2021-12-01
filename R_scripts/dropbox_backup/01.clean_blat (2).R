### INFO: 
### DATE: Mon Jul 08 21:52:06 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/flanking_genes.mRNA/results/blat_Elob_psgene")

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

# blat main path
blat_path <- file.path(inpath)

# blat .psl path
blat_psl_path <- list.files(path = blat_path, pattern = ".*\\.score\\.psl$", full.names = T)

# blat .bed path (with strand info)
blat_bed_path <- list.files(path = blat_path, pattern = ".*\\.bed$", full.names = T)

# fasta index path
fasta_index_path <- list.files(path = blat_path, pattern = ".*\\.fa\\.fai$", full.names = T)

# inter exons bed 
inter_exons_bed_path <- list.files(path = file.path(inpath, ".."), pattern = ".*\\.inter_exons\\.bed", full.names = T)

######################################################## READ DATA
# read blat .psl
blat_psl <- 
  purrr::map(blat_psl_path, function(path){
    
    # read and set column names of .psl, arrange by score
    suppressMessages(readr::read_delim(file = path, delim = "\t", col_names = c("seqnames", "start", "end", "subject_coords", "score", "perc_identity")) %>% 
                       dplyr::mutate(start = start + 1, 
                                     species = basename(path) %>% str_remove(., "\\.score\\.psl"), 
                                     query_coords = str_c(species, ":", seqnames, ":", start, "-", end)))
    
  }) %>% 
  dplyr::bind_rows(.)

# read blat .bed
blat_bed <- 
  purrr::map(blat_bed_path, function(path){
    
    # read bed
    suppressMessages(readr::read_delim(path, delim = "\t", 
                                       col_names = c("seqnames", "start", "end", "name", "score", "strand", 
                                                     "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", 
                                                     "blockStarts")) %>% 
                       dplyr::select(seqnames, start, end, strand) %>% 
                       dplyr::mutate(start = start + 1, 
                                     width = (end - start) + 1, 
                                     species = basename(path) %>% str_remove(., "\\.bed"), 
                                     query_coords = str_c(species, ":", seqnames, ":", start, "-", end)) %>% 
                       dplyr::select(query_coords, width, strand))
    
  }) %>% 
  dplyr::bind_rows(.)

# read fasta index
fasta_index <- 
  purrr::map(fasta_index_path, function(path){
    
    # read index, add species
    suppressMessages(readr::read_delim(file = path, delim = "\t", col_name = c("transcript_id", "mm10_total_width", "tmp1", "tmp2", "tmp3")) %>% 
                       dplyr::mutate(gene_name = basename(path) %>% str_remove(., "\\.fa\\.fai")) %>% 
                       dplyr::select(gene_name, mm10_total_width))
    
  }) %>% 
  dplyr::bind_rows(.)

# read inter exons .bed
inter_exons_bed <- 
  purrr::map(inter_exons_bed_path, function(path){
    
    # read bed
    suppressMessages(readr::read_delim(path, delim = "\t", 
                                       col_names = c("seqnames", "start", "end", "name", "score", "strand", 
                                                     "thickStart", "thickEnd", "itemRgb")) %>% 
                       dplyr::select(seqnames, start, end, strand) %>% 
                       dplyr::mutate(start = start + 1, 
                                     width = (end - start) + 1, 
                                     species = basename(path) %>% str_remove(., "\\.inter_exons\\.bed"), 
                                     query_coords = str_c(species, ":", seqnames, ":", start, "-", end)))
  }) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# add strand info, get top score for each species
blat_results_tb <- 
  blat_psl %>% 
  dplyr::left_join(., blat_bed, by = "query_coords") %>%
  dplyr::mutate(subject_coords = str_remove(subject_coords, ".*:")) %>% 
  tidyr::separate(subject_coords, c("hit_start", "hit_end"), sep = "-") %>%
  dplyr::mutate(hit_start = as.numeric(hit_start) + 1, 
                hit_end = as.numeric(hit_end),
                mm10_hit_width = (hit_end - hit_start) + 1) %>% 
  dplyr::select(seqnames, start, end, width, strand, score, perc_identity, mm10_hit_width, species) %>% 
  dplyr::group_by(species) %>% 
  top_n(1, score) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(gene_name = str_extract(species, str_c(fasta_index$gene_name, collapse = "|")), 
                species = str_remove(species, str_c("\\.", gene_name))) %>% 
  dplyr::left_join(., fasta_index, by = "gene_name") %>% 
  dplyr::mutate(mm10_hit_perc = round((mm10_hit_width / mm10_total_width), 3) * 100) %>% 
  dplyr::select(species, gene_name, everything()) 

y <- 
  blat_results_tb %>% 
  right_join(., inter_exons_bed, by = c("species", "seqnames")) %>% 
  dplyr::select(species, seqnames, everything())

# # split to list of species
# blat_results_filt <- 
#   blat_results_tb %>% 
#   split(., .$species)
# 
# # save as .bed
# purrr::map(names(blat_results_filt), function(species){
#   
#   # filter by species
#   tb <- blat_results_filt[[species]]
#   
#   # filter by gene name, save as .bed
#   purrr::map(unique(tb$gene_name), function(gene){
#     
#     # save as .bed
#     tb %>% 
#       dplyr::filter(gene_name == gene) %>% 
#       GRanges(.) %T>%
#       rtracklayer::export.bed(., file.path(outpath, "results", str_c(species, gene, "top.bed", sep = ".")))
#     
#   })
#   
#   
# })


# # save as table
# bind_rows(blat_results_filt) %>% 
#   tidyr::unite(coords, seqnames, start, end, sep = " ") %T>%
#   readr::write_csv(., file.path(outpath, "mm10.lnc1_3prime.top_results.blat.csv"))


