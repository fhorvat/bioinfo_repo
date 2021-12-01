### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/final_cluster_expression")

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

# RPM table path
rpm_tb_path <- file.path(inpath, "PIWIL_IP")
rpm_tb_path <- list.files(rpm_tb_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- list.files(inpath, pattern = "piRNA_clusters.*\\.norm_width\\.csv", full.names = T)

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

######################################################## MAIN CODE
# set table name
table_name <- 
  rpm_tb_path %>% 
  basename(.) 

# get clusters and their RPM values as GenomicRanges
rpm_gr <- 
  rpm_tb %>% 
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  dplyr::rename(norm_width = width)

# get clusters RPKM
clusters_rpkm <- 
  rpm_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(coordinates, PIWIL1, PIWIL3, norm_width) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N) 

# save 
readr::write_csv(clusters_rpkm, file.path(outpath, table_name %>% str_replace(., "FPM_mean.csv", "RPKM.csv")))
