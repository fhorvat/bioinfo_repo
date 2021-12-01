### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path <- list.files(expression_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# normalized width path
norm_width_path <- list.files(inpath, pattern = "\\.norm_width\\.csv", full.names = T, recursive = T)

# clusters path
clusters_path <- list.files(inpath, ".*\\.merged_clusters\\.csv", full.names = T)

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

# read cluster coordinates
clusters_tb_list <- purrr::map(clusters_path, function(path){
  
  # read table, keep only coordinates, separate
  clusters_tmp <- 
    readr::read_csv(path) %>% 
    dplyr::select(coordinates) %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") 
  
  # return 
  return(clusters_tmp)
  
}) %>% 
  set_names(., str_extract(clusters_path, "prepachytene|pachytene"))

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.1k_windows"

# set RPKM cutoff
rpkm_cutoff <- 1

# get original tiles as GenomicRanges
rpm_gr <- 
  rpm_tb %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>%
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  dplyr::rename(norm_width = width) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# overlap with ranges
clusters_rpm_list <- purrr::map(names(clusters_tb_list), function(stage){
  
  # get table, convert to genomicRanges
  clusters_gr <- 
    clusters_tb_list[[stage]] %>% 
    GRanges(.)
  
  # overlap
  overlaps <- GenomicRanges::findOverlaps(rpm_gr, clusters_gr, ignore.strand = T)
  
  # get unique clusters table
  clusters_tb <- 
    clusters_gr[subjectHits(overlaps)] %>% 
    as_tibble(.) %>% 
    dplyr::select(-c(strand, width)) %>% 
    dplyr::mutate(cluster_id = subjectHits(overlaps)) %>% 
    unique(.)
  
  # get clusters total RPKM
  clusters_rpkm <- 
    rpm_gr[queryHits(overlaps)] %>% 
    as_tibble(.) %>% 
    dplyr::mutate(cluster_id = subjectHits(overlaps)) %>% 
    dplyr::select(cluster_id, WT_13dpp, KO_13dpp, WT_21dpp, KO_21dpp, norm_width) %>% 
    dplyr::group_by(cluster_id) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::left_join(., clusters_tb, by = "cluster_id") %>% 
    dplyr::select(-cluster_id) %>% 
    tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
    tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
    dplyr::mutate(norm_width = (norm_width / 1000), 
                  rpkm = round((rpm / norm_width), 3)) %>% 
    tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = "rpkm") %>% 
    dplyr::mutate(log2FC_13dpp = log2(KO_13dpp / WT_13dpp), 
                  log2FC_21dpp = log2(KO_21dpp / WT_21dpp)) %>% 
    dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0)))
  
  # return
  return(clusters_rpkm)
  
}) %>% 
  set_names(., names(clusters_tb_list))


## filter and save
# prepachytene
prepachytene_tb <- 
  clusters_rpm_list[["prepachytene"]] %>% 
  dplyr::filter(log2FC_13dpp <= -2, log2FC_21dpp <= -2) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(width = width / 1000) %>% 
  dplyr::rename(width_kb = width, norm_width_kb = norm_width) %>% 
  dplyr::arrange(-norm_width_kb) %T>% 
  readr::write_csv(., file.path(outpath, "MesAur1.1k_windows.rpkm_cutoff.1.preachytene.merged_clusters.RPKM.csv"))

# pachytene
pachytene_tb <- 
  clusters_rpm_list[["pachytene"]] %>% 
  dplyr::filter(log2FC_13dpp > -2, log2FC_21dpp <= -2) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(width = width / 1000) %>% 
  dplyr::rename(width_kb = width, norm_width_kb = norm_width) %>% 
  dplyr::arrange(-norm_width_kb) %T>% 
  readr::write_csv(., file.path(outpath, "MesAur1.1k_windows.rpkm_cutoff.1.pachytene.merged_clusters.RPKM.csv"))
