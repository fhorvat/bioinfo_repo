### INFO:
### DATE: Mon Jul 08 21:52:06 2019
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

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"
  
# blat .psl path
blat_psl_path <- list.files(path = inpath, pattern = ".*\\.score\\.psl$", full.names = T)

# blat .bed path (with strand info)
blat_bed_path <- list.files(path = inpath, pattern = ".*\\.bed$", full.names = T)

# fasta index path
fasta_index_path <- list.files(path = genome_path, pattern = ".*\\.fasta\\.fai$", full.names = T)

######################################################## READ DATA
# read blat .psl
blat_psl <- 
  readr::read_delim(file = blat_psl_path, 
                    delim = "\t", 
                    col_names = c("seqnames", "start", "end", "subject_coords", "score", "perc_identity")) %>%
  dplyr::mutate(start = start + 1,
                query_coords = str_c(seqnames, ":", start, "-", end))

# read blat .bed
blat_bed <- 
  readr::read_delim(file = blat_bed_path, 
                    delim = "\t",
                    col_names = c("seqnames", "start", "end", "name", "score", "strand",
                                  "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                                  "blockStarts")) %>%
  dplyr::select(seqnames, start, end, strand) %>%
  dplyr::mutate(start = start + 1,
                width = (end - start) + 1,
                query_coords = str_c(seqnames, ":", start, "-", end)) %>%
  dplyr::select(query_coords, width, strand)

# read fasta index
fasta_index <- 
  readr::read_delim(file = fasta_index_path, delim = "\t", col_name = c("transcript_id", "genome_total_width", "tmp1", "tmp2", "tmp3")) %>%
  dplyr::select(transcript_id, genome_total_width)

######################################################## MAIN CODE
# add strand info, get top score for each species
blat_results_tb <-
  blat_psl %>%
  dplyr::left_join(., blat_bed, by = "query_coords") %>%
  dplyr::mutate(subject_coords = str_remove(subject_coords, ".*:")) %>%
  tidyr::separate(subject_coords, c("hit_start", "hit_end"), sep = "-") %>%
  dplyr::mutate(hit_start = as.numeric(hit_start) + 1,
                hit_end = as.numeric(hit_end),
                genome_hit_width = (hit_end - hit_start) + 1) %>%
  dplyr::select(seqnames, start, end, width, strand, score, perc_identity, genome_hit_width) %>% 
  dplyr::arrange(-score)
