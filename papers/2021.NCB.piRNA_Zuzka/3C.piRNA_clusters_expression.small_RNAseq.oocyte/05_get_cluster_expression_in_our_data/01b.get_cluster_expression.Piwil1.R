### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/oocyte_expression/piwil1")

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
rpm_tb_path <- file.path(inpath, "feature_counts")
rpm_tb_path <- list.files(rpm_tb_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- "../../count_N_in_genome"
norm_width_path <- list.files(norm_width_path, 
                              pattern = basename(rpm_tb_path) %>% str_replace(., "\\.FPM_mean\\.csv$", ".norm_width.csv"), 
                              full.names = T)

# cluster data path
cluster_path <- file.path(inpath, 
                          "../../recalculate_cluster_expression",
                          "piRNA_clusters.oocyte_PIWIL1.rpkm_cutoff.1.20210513.RPKM.csv")

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

# read cluster expression table
cluster_tb <- readr::read_csv(cluster_path)

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
  dplyr::select(coordinates, starts_with("Mov10l1_"), norm_width) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
  # dplyr::mutate(rpkm_log2FC_PIWIL1_PIWIL3 = log2(rpkm_PIWIL1 / rpkm_PIWIL3), 
  #               rpm_log2FC_PIWIL1_PIWIL3 = log2(rpm_PIWIL1 / rpm_PIWIL3)) %>%  
  # dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N) 

# join with original cluster data
clusters_rpkm %>% 
  dplyr::left_join(., cluster_tb %>% select(-c(width_full, width_without_N)), by = "coordinates") %>% 
  dplyr::select(coordinates:width_without_N, contains("PIWIL"), everything()) %>% 
  readr::write_csv(., file.path(outpath, table_name %>% str_replace(., "FPM_mean.csv", "RPKM.deduplexed.our_data.csv")))




