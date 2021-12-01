### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_homologs")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# blast results path
blast_results_path <- list.files(inpath, "Mov10.protein.tblastn.txt", full.names = T)

######################################################## READ DATA
# read blast results
blast_results <- readr::read_delim(file = blast_results_path, 
                                   delim = "\t", 
                                   col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                 "mismatches", "gap_open", 
                                                 "query_start", "query_end", "subject_start", "subject_end", 
                                                 "e_value", "bit_score")) 

######################################################## MAIN CODE
# filter results
blast_hits <- 
  blast_results %>% 
  dplyr::mutate(strand = ifelse(subject_start < subject_end, "+", "-")) %>% 
  dplyr::mutate(start = ifelse(strand == "+", subject_start, subject_end), 
                end = ifelse(strand == "+", subject_end, subject_start)) %>% 
  dplyr::select(seqnames = subject_id, start, end, strand, query_start, query_end, identity_perc, bit_score) %>% 
  dplyr::mutate(query_alignment_width = query_end - query_start + 1, 
                subject_alignment_width = end - start + 1)

# to GRanges
blast_gr <- GRanges(blast_hits)

# reduce close hits (exons) to genes
blast_gr_reduced <- GenomicRanges::reduce(blast_gr, ignore.strand = F, min.gapwidth = 20000)
mcols(blast_gr_reduced)$gene_id <- str_c(seqnames(blast_gr_reduced), start(blast_gr_reduced), end(blast_gr_reduced), sep = " ")

# for each exon get cooresponding gene
overlaps <- findOverlaps(blast_gr, blast_gr_reduced)

# add gene info to exons
exons_tb <- 
  blast_gr[queryHits(overlaps)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(gene_hit =  blast_gr_reduced[subjectHits(overlaps)] %>% mcols %$% gene_id) %>% 
  dplyr::group_by(gene_hit) %>% 
  dplyr::summarise(subject_alignment_width = sum(subject_alignment_width), 
                   bit_score = sum(bit_score)) %>% 
  dplyr::arrange(-bit_score)
