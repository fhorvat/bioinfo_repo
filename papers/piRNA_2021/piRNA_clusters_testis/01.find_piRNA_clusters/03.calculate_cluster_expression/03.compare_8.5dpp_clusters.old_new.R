### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/rpkm_cutoff_10.clusters.Petr.20210104")

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
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters path
old_path <- file.path(inpath, "Table S1 9 dpp pre-pachytene piRNA clusters.xlsx")
new_path <- file.path(inpath, "Table S1 9 dpp pre-pachytene piRNA clusters 210104.xlsx")

######################################################## READ DATA
# read clusters
clusters_old <- openxlsx::read.xlsx(old_path) %>% as_tibble(.)
clusters_new <- openxlsx::read.xlsx(new_path) %>% as_tibble(.)

######################################################## MAIN CODE
# get coordinates - old clusters 
clusters_old_gr <- 
  clusters_old %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# get coordinates - new clusters
clusters_new_gr <- 
  clusters_new %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# overlap
overlaps <- findOverlaps(clusters_new_gr, clusters_old_gr)

# get all clusters overlaping with old clusters
clusters_new_gr[unique(queryHits(overlaps))]

length(clusters_new_gr)
length(unique(queryHits(overlaps)))

length(clusters_old_gr)
length(unique(subjectHits(overlaps)))


