### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/count_N_in_genome")

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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window.no_multimappers")

# RPM table path
rpm_tb_path <- list.files(expression_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# get normalized tile width path
tile_width_path <- list.files(inpath, "MesAur1.1k_windows.norm_width.csv", full.names = T)

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized tile width table
tile_tb <- readr::read_csv(tile_width_path)

######################################################## MAIN CODE
# normalize rpm
rpkm_tb <-
  rpm_tb %>%
  tidyr::pivot_longer(-c(gene_id, coordinates), values_to = "rpm", names_to = "sample_id") %>%
  dplyr::left_join(., tile_tb, by = "gene_id") %>%
  dplyr::filter(width > 0) %>%
  dplyr::mutate(width = (width / 1000),
                rpkm = round((rpm / width), 3)) %>%
  tidyr::pivot_wider(id_cols = c(gene_id, coordinates), names_from = "sample_id", values_from = "rpkm")

# save
readr::write_csv(rpkm_tb, file.path(inpath, rpm_tb_path %>% basename(.) %>% str_replace(., "FPM", "8.5dpp.no_multimappers.RPKM")))


