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
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.manual_clusters")

# RPM table path
rpm_tb_path <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- list.files(inpath, pattern = "\\.norm_width\\.csv", full.names = T, recursive = T)

# 13dpp and 21dpp results path
results_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.csv")
results_path.reseq <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.reseq.csv")

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_tb_path) %>% 
  as_tibble(.)

# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

# read 13dpp and 21dpp results 
results_tb <- readr::read_csv(results_path)
results_tb.reseq <- readr::read_csv(results_path.reseq)

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
  dplyr::select(coordinates, WT_8.5dpp = WT_8.5, KO_8.5dpp = KO_8.5, norm_width) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
  dplyr::mutate(rpkm_log2FC_8.5dpp = log2(rpkm_KO_8.5dpp / rpkm_WT_8.5dpp), 
                rpm_log2FC_8.5dpp = log2(rpm_KO_8.5dpp / rpm_WT_8.5dpp)) %>% 
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
readr::write_csv(clusters_tb_sorted, file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.8.5dpp.csv"))

#### join with 13dpp and 21dpp results and save
# rename reseq table
results_tb.reseq %<>% 
  dplyr::select(-c(width_full, width_without_N)) %>% 
  dplyr::rename_at(vars(starts_with("rp")), ~str_c(., ".reseq"))

# join with all results and save
clusters_tb_sorted %>% 
  dplyr::left_join(., results_tb %>% dplyr::select(-c(width_full, width_without_N)), by = "coordinates") %>% 
  dplyr::left_join(., results_tb.reseq, by = "coordinates") %T>%
  readr::write_csv(., file.path(outpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.8.5dpp.13dpp.21dpp.csv"))

