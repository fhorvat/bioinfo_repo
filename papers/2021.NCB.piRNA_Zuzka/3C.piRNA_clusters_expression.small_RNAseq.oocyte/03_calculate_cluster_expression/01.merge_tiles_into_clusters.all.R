### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/small_RNAseq")

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
### get clusters
getClusters <- function(tiles_tb, super_clusters = T, merge_width = 2000){
  
  # get coordinates, merge neighbours
  clusters_gr <- 
    tiles_tb %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
    GRanges(.) %>% 
    GenomicRanges::reduce(., ignore.strand = T)
  
  # merge to super clusters
  if(super_clusters){
    
    # join clusters which are up to 2kb appart
    clusters_gr <- GenomicRanges::reduce(clusters_gr, min.gapwidth = merge_width)
    
  }
  
  # overlap with original tiles
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
    dplyr::select(cluster_id, oocyte_PIWIL1, oocyte_PIWIL3, norm_width) %>%
    dplyr::group_by(cluster_id) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::left_join(., clusters_tb, by = "cluster_id") %>% 
    dplyr::select(-cluster_id) %>% 
    tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
    tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
    dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
    tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
    dplyr::mutate(rpkm_log2FC.PIWIL1_PIWIL3 = log2(rpkm_oocyte_PIWIL1 / rpkm_oocyte_PIWIL3), 
                  rpm_log2FC.PIWIL1_PIWIL3 = log2(rpm_oocyte_PIWIL1 / rpm_oocyte_PIWIL3)) %>% 
    dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0)))
  
  # return table
  return(clusters_rpkm)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### Piwil1/Piwil3 IP
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_Siomi_WT.Piwil1_Piwil3_IP.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# RPKM table path
rpkm_tb_path <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.RPKM_mean.Piwil1_Piwil3_IP.small_RNAseq.csv")


### normalized widths
# normalized width path
norm_width_path <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.norm_width.csv")

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path) 

# read RPKM table
rpkm_tb <- readr::read_csv(rpkm_tb_path)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.1k_windows"

# set RPKM cutoff
rpkm_cutoff <- 1

# get original tiles and their RPM values as GenomicRanges
rpm_gr <- 
  rpm_tb %>% 
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  dplyr::rename(norm_width = width) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)


### get clusters 
# calculate fold change in 1 kb tiles
tiles_fc <-
  rpkm_tb %>%
  # set_colnames(., str_remove(colnames(.), "Mov10l1_")) %>%
  # dplyr::filter_at(.vars = vars(contains("P")), .vars_predicate = any_vars(. > rpkm_cutoff)) %>%
  dplyr::mutate(log2FC.piwil1_piwil3 = log2(oocyte_PIWIL1 / oocyte_PIWIL3))


## PIWIL1 clusters 
# merge those who pass RPKM cutoff, merge again if clusters are at most 2kb apart
clusters_piwil1 <- 
  tiles_fc %>% 
  dplyr::filter(oocyte_PIWIL1 > rpkm_cutoff) %>%
  getClusters(tiles_tb = ., super_clusters = T, merge_width = 2001) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N)

clusters_piwil1 %>% 
  dplyr::select(coordinates, width_full, width_without_N, rpkm_oocyte_PIWIL1, rpm_oocyte_PIWIL1) %T>% 
  readr::write_csv(., file.path(outpath, "piRNA_clusters.oocyte_PIWIL1.rpkm_cutoff.1.20210513.csv"))

## PIWIL3 clusters 
# merge those who pass RPKM cutoff, merge again if clusters are at most 2kb apart
clusters_piwil3 <- 
  tiles_fc %>% 
  dplyr::filter(oocyte_PIWIL3 > rpkm_cutoff) %>%
  getClusters(tiles_tb = ., super_clusters = T, merge_width = 2001) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N)

clusters_piwil3 %>% 
  dplyr::select(coordinates, width_full, width_without_N, rpkm_oocyte_PIWIL3, rpm_oocyte_PIWIL3) %T>% 
  readr::write_csv(., file.path(outpath, "piRNA_clusters.oocyte_PIWIL3.rpkm_cutoff.1.20210513.csv"))

