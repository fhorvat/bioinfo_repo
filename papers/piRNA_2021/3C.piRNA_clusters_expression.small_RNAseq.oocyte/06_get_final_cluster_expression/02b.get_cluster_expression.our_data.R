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
rpm_tb_path <- file.path(inpath, "in_house/feature_counts")
rpm_tb_path <- list.files(rpm_tb_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- file.path(inpath, "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.norm_width.csv")

# cluster expression path
cluster_path <- file.path(inpath, "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.RPKM.csv")

# cluster annotation path
cluster_annotation_path <- file.path(inpath, "../oocyte_expression/piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.csv")

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

# read cluster expression table
cluster_tb <- readr::read_csv(cluster_path)

# read cluster annotation path
cluster_annotation <- readr::read_csv(cluster_annotation_path)

######################################################## MAIN CODE
# tidy cluster annotation
cluster_annotation %<>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(coordinates, cluster_in)

# set table name
table_name <- "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210524.final.csv"

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
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N) 

### join with original cluster data
# coordinates MesAur 1.0	width full	width without Ns	
# rpkm oocyte	 rpm oocyte	 
# piRNA type	
# rpkm PIWIL1 IP	rpkm PIWIL3 IP	rpm PIWIL1 IP	rpm PIWIL3 IP
clusters_rpkm_final <- 
  clusters_rpkm %>% 
  dplyr::left_join(., cluster_tb %>% select(-c(width_full, width_without_N)), by = "coordinates") %>% 
  dplyr::left_join(., cluster_annotation, by = "coordinates") %>% 
  dplyr::mutate(rpkm_oocyte = rpkm_Mov10l1_WT.18to20nt + rpkm_Mov10l1_WT.24to31nt, 
                rpm_oocyte = rpm_Mov10l1_WT.18to20nt + rpm_Mov10l1_WT.24to31nt) %>% 
  dplyr::select(`coordinates MesAur 1.0` = coordinates, `width full` = width_full, `width without Ns` = width_without_N, 
                `rpkm oocyte` = rpkm_oocyte, `rpm oocyte` = rpm_oocyte, 
                `piRNA type` = cluster_in, 
                `rpkm PIWIL1 IP` = rpkm_PIWIL1, `rpkm PIWIL3 IP` = rpkm_PIWIL3, 
                `rpm PIWIL1 IP` = rpm_PIWIL1, `rpm PIWIL3 IP` = rpm_PIWIL3) %>% 
  dplyr::filter_at(vars(starts_with("rpm PIWIL")), dplyr::any_vars(. >= 10)) %>% 
  dplyr::filter(`rpm oocyte` >= 10) %>% 
  dplyr::arrange(-`width without Ns`) %T>% 
  readr::write_csv(., file.path(outpath, table_name))




