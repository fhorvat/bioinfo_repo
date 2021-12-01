### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/recalculate_clusters_RPM_and_RPKM")

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

library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# input table path
clusters_tb_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters_final_full_reduced_200730-for recalculation.xlsx")

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_tb_path) %>% 
  as_tibble(.)

######################################################## MAIN CODE
# get coordinates
clusters_gr <-
  clusters_tb %>% 
  dplyr::select(coordinates) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# set names
names(clusters_gr) <- str_c(seqnames(clusters_gr), ":", start(clusters_gr), "-", end(clusters_gr))

# export as .gtf
export.gff3(object = clusters_gr, con = file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.gff"))

# export as .bed
export.bed(object = clusters_gr, con = file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.bed"))


# # save as SAF
# line1_saf <-
#   line1_gr %>%
#   as_tibble(.) %>%
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>%
#   readr::write_delim(., file.path(outpath, "MesAur1.5k_windows.saf"), delim = "\t")

