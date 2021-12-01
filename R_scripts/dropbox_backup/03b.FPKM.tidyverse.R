### INFO: 
### DATE: Wed Oct 30 07:36:10 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# counts table path
counts_path <- file.path(inpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.chr5.csv")

# sample table path
sample_table_path <- file.path(inpath, "CNOT6L.sample_table.csv")
  
# genes info path
genes_info_path <- file.path(inpath, "ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.chr5_chr6.csv")

######################################################## READ DATA
# read tables
counts_tb <- readr::read_csv("ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.chr5.csv")
sample_table <- readr::read_csv("CNOT6L.sample_table.csv")
genes_info_tb <- readr::read_csv("ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.chr5_chr6.csv")

######################################################## MAIN CODE
### FPKM
# get long table of counts, join with tables, calculate FPKMs, return to wide
fpkm_tb <-
  counts_tb %>% 
  tidyr::pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "counts") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.Aligned\\.sortedByCoord\\.out\\.bam$")) %>% 
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., genes_info_tb %>% dplyr::select(gene_id, width = total_exon_length_sum), by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  tidyr::pivot_wider(names_from = sample_id, values_from = fpkm)
