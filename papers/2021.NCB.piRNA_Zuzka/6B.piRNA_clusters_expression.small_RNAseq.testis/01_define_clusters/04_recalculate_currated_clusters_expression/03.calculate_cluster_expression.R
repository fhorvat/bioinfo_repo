### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
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

library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# input table path
clusters_tb_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters_final_full_reduced_200730-for recalculation.xlsx")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.manual_clusters")

# RPM table path
rpm_tb_path <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- list.files(inpath, pattern = "\\.norm_width\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_tb_path) %>% 
  as_tibble(.)

# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.1k_pachytene_clusters.200730"

# get clusters and their RPM values as GenomicRanges
rpm_gr <- 
  rpm_tb %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>%
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  dplyr::rename(norm_width = width)

# get clusters RPKM
clusters_rpkm <- 
  rpm_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(coordinates, WT_13dpp, KO_13dpp, WT_21dpp, KO_21dpp, norm_width) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
  dplyr::mutate(rpkm_log2FC_13dpp = log2(rpkm_KO_13dpp / rpkm_WT_13dpp), 
                rpkm_log2FC_21dpp = log2(rpkm_KO_21dpp / rpkm_WT_21dpp), 
                rpm_log2FC_13dpp = log2(rpm_KO_13dpp / rpm_WT_13dpp), 
                rpm_log2FC_21dpp = log2(rpm_KO_21dpp / rpm_WT_21dpp)) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N) 

# arrange as original table
clusters_tb_sorted <- 
  clusters_tb %>% 
  dplyr::select(coordinates) %>% 
  left_join(., clusters_rpkm, by = "coordinates") 

# save 
readr::write_csv(clusters_tb_sorted, file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.csv"))
