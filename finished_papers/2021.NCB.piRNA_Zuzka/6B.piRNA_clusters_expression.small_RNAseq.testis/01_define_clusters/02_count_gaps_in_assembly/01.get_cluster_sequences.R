### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/count_N_in_genome")

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
library(Biostrings)
library(BSgenome.Maur.UCSC.MesAur1)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# .gff file path
gff_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Analysis/expression.1k_window"
gff_path <- file.path(gff_path, "MesAur1.1k_windows.gff")

######################################################## READ DATA
# read gff
gff_tb <- readr::read_delim(gff_path, delim = "\t", comment = "#", col_names = F)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.1k_windows"

# get GRanges coordinates
tile_gr <- 
  gff_tb %>% 
  dplyr::select(seqnames = X1, start = X4, end = X5) %>% 
  GRanges(.)
names(tile_gr) <- str_c(seqnames(tile_gr), ":", start(tile_gr), "-", end(tile_gr))

# get tile sequences 
tile_seq <- getSeq(BSgenome.Maur.UCSC.MesAur1, tile_gr)

# save
Biostrings::writeXStringSet(tile_seq, file.path(outpath, str_c(table_name, "fa", sep = ".")))


