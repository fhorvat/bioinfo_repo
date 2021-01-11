### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/map_clusters_to_Siomi")

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
library(BSgenome.Maur.UCSC.MesAur1)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters table path
clusters_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.xlsx")

######################################################## READ DATA
# read cluster tables
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(., .name_repair = "unique")

######################################################## MAIN CODE
# get coordinates
clusters_gr <- 
  clusters_tb %>% 
  dplyr::select(coordinates) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
  dplyr::filter(!is.na(seqnames)) %>% 
  GRanges(.)

# get sequences
clusters_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.MesAur1, clusters_gr)
names(clusters_seq) <- str_replace_all(mcols(clusters_gr)$coordinates, " ", "_")

# save
Biostrings::writeXStringSet(clusters_seq, file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.fasta"))


