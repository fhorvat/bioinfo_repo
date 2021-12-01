### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/map_clusters_to_Siomi")

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
clusters_path <- file.path(inpath, "../expression/piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.csv")

######################################################## READ DATA
# read cluster tables
clusters_tb <- readr::read_csv(clusters_path)

######################################################## MAIN CODE
# get coordinates
clusters_gr <- 
  clusters_tb %>% 
  dplyr::select(coordinates, cluster_id) %>% 
  dplyr::mutate(cluster_id = str_replace_all(cluster_id, " ", "_")) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = T) %>% 
  dplyr::filter(!is.na(seqnames)) %>% 
  GRanges(.)

# get sequences
clusters_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.MesAur1, clusters_gr)
names(clusters_seq) <- mcols(clusters_gr)$cluster_id

# save
Biostrings::writeXStringSet(clusters_seq, file.path(outpath, "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.fasta"))


