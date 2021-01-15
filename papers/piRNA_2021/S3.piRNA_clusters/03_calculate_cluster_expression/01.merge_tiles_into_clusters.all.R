### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq")

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
    dplyr::select(cluster_id, WT.21dpp:KO.8.5dpp, norm_width) %>%
    dplyr::group_by(cluster_id) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::left_join(., clusters_tb, by = "cluster_id") %>% 
    dplyr::select(-cluster_id) %>% 
    tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
    tidyr::pivot_longer(cols = -c(coordinates, norm_width), values_to = "rpm", names_to = "sample_id") %>% 
    dplyr::mutate(rpkm = round((rpm / (norm_width / 1000)), 3)) %>% 
    tidyr::pivot_wider(id_cols = c(coordinates, norm_width), names_from = "sample_id", values_from = c("rpkm", "rpm")) %>% 
    dplyr::mutate(rpkm_log2FC.8.5dpp = log2(rpkm_KO.8.5dpp / rpkm_WT.8.5dpp), 
                  rpm_log2FC.8.5dpp = log2(rpm_KO.8.5dpp / rpm_WT.8.5dpp), 
                  rpkm_log2FC.13dpp = log2(rpkm_KO.13dpp / rpkm_WT.13dpp), 
                  rpm_log2FC.13dpp = log2(rpm_KO.13dpp / rpm_WT.13dpp), 
                  rpkm_log2FC.21dpp = log2(rpkm_KO.21dpp / rpkm_WT.21dpp), 
                  rpm_log2FC.21dpp = log2(rpm_KO.21dpp / rpm_WT.21dpp)) %>% 
    dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0)))
  
  # return table
  return(clusters_rpkm)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### 21 dpp path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path_21dpp <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# RPKM table path
rpkm_tb_path_21dpp <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.RPKM_mean.csv")


### 13 dpp path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path_13dpp <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# RPKM table path
rpkm_tb_path_13dpp <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.RPKM_mean.13dpp.reseq.csv")


### 8.5 dpp path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path_8.5dpp <- list.files(expression_path, pattern = ".*\\.FPM_mean\\.csv$", full.names = T)

# RPKM table path
rpkm_tb_path_8.5dpp <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.RPKM_mean.8.5dpp.run_2.csv")


### normalized widths
# normalized width path
norm_width_path <- file.path(inpath, "../count_N_in_genome/MesAur1.1k_windows.norm_width.csv")

######################################################## READ DATA
### 21 dpp
# read RPM table
rpm_tb_21dpp <- 
  readr::read_csv(rpm_tb_path_21dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.21dpp = Mov10l_WT_21dpp, Mov10l1_KO.21dpp = Mov10l_KO_21dpp)

# read RPKM table
rpkm_tb_21dpp <- 
  readr::read_csv(rpkm_tb_path_21dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.21dpp = Mov10l_WT_21dpp, Mov10l1_KO.21dpp = Mov10l_KO_21dpp)


### 13 dpp
# read RPM table
rpm_tb_13dpp <- 
  readr::read_csv(rpm_tb_path_13dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.13dpp = Mov10l_WT_13dpp, Mov10l1_KO.13dpp = Mov10l_KO_13dpp)

# read RPKM table
rpkm_tb_13dpp <- 
  readr::read_csv(rpkm_tb_path_13dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.13dpp = Mov10l_WT_13dpp, Mov10l1_KO.13dpp = Mov10l_KO_13dpp)


### 8.5 dpp
# read RPM table
rpm_tb_8.5dpp <- 
  readr::read_csv(rpm_tb_path_8.5dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.8.5dpp = Mov10l1_WT_8.5dpp, Mov10l1_KO.8.5dpp = Mov10l1_KO_8.5dpp)

# read RPKM table
rpkm_tb_8.5dpp <- 
  readr::read_csv(rpkm_tb_path_8.5dpp) %>% 
  dplyr::select(gene_id, coordinates, Mov10l1_WT.8.5dpp = Mov10l1_WT_8.5dpp, Mov10l1_KO.8.5dpp = Mov10l1_KO_8.5dpp)


### other data
# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

######################################################## MAIN CODE
### clean and join data
# RPM
rpm_tb <- 
  rpm_tb_21dpp %>% 
  dplyr::left_join(., rpm_tb_13dpp, by = c("gene_id", "coordinates")) %>% 
  dplyr::left_join(., rpm_tb_8.5dpp, by = c("gene_id", "coordinates"))

# RPKM
rpkm_tb <- 
  rpkm_tb_21dpp %>% 
  dplyr::left_join(., rpkm_tb_13dpp, by = c("gene_id", "coordinates")) %>% 
  dplyr::left_join(., rpkm_tb_8.5dpp, by = c("gene_id", "coordinates"))

# set name
table_name <- "MesAur1.1k_windows"

# set RPKM cutoff
rpkm_cutoff <- 1

# get original tiles and their RPM values as GenomicRanges
rpm_gr <- 
  rpm_tb %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l1_")) %>%
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  dplyr::rename(norm_width = width) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)


### get clusters 
# calculate fold change in 1 kb tiles
tiles_fc <-
  rpkm_tb %>%
  set_colnames(., str_remove(colnames(.), "Mov10l1_")) %>%
  dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > rpkm_cutoff)) %>%
  dplyr::mutate(log2FC.8.5dpp = log2(KO.8.5dpp / WT.8.5dpp), 
                log2FC.13dpp = log2(KO.13dpp / WT.13dpp), 
                log2FC.21dpp = log2(KO.21dpp / WT.21dpp))

# 8.5 dpp clusters = merge if log2FC 8.5D KO/WT < -2
# merge again if clusters are at most 2kb apart
clusters_8.5dpp <- 
  tiles_fc %>% 
  dplyr::filter(log2FC.8.5dpp <= -2) %>% 
  getClusters(tiles_tb = ., super_clusters = T, merge_width = 2001) %>% 
  dplyr::filter(rpkm_log2FC.8.5dpp <= -2) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N)
clusters_8.5dpp %>% 
  dplyr::select(coordinates, width_full, width_without_N, rpkm_8.5dpp = rpkm_WT.8.5dpp, rpm_8.5dpp = rpm_WT.8.5dpp) %T>% 
  readr::write_csv(., file.path(outpath, "piRNA_clusters.8.5dpp.rpkm_cutoff.1.20201231.csv"))

# prepachytene clusters = merge if log2FC 13D KO/WT < -2 and log2FC 21D < -2
# merge again if clusters are at most 2kb apart
clusters_13dpp <- 
  tiles_fc %>% 
  dplyr::filter(log2FC.13dpp <= -2, log2FC.21dpp <= -2) %>% 
  getClusters(tiles_tb = ., super_clusters = T, merge_width = 2001) %>% 
  dplyr::filter(rpkm_log2FC.13dpp <= -2, rpkm_log2FC.21dpp <= -2) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>% 
  dplyr::arrange(-width_without_N) 
clusters_13dpp %>% 
  dplyr::select(coordinates, width_full, width_without_N, rpkm_13dpp = rpkm_WT.13dpp, rpm_13dpp = rpm_WT.13dpp) %T>% 
  readr::write_csv(., file.path(outpath, "piRNA_clusters.13dpp.rpkm_cutoff.1.20201231.csv"))

# pachytene clusters = merge if log2FC 21D KO/WT <-2, log2FC 13D KO/WT >= -2
# merge again if clusters are at most 1kb appart
clusters_21dpp <- 
  tiles_fc %>% 
  dplyr::filter(log2FC.13dpp > -2, log2FC.21dpp <= -2) %>% 
  getClusters(tiles_tb = ., super_clusters = T, merge_width = 1001) %>% 
  dplyr::filter(rpkm_log2FC.13dpp > -2, rpkm_log2FC.21dpp <= -2) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>%  
  as_tibble(.) %>% 
  select(-strand) %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::rename(width_full = width, width_without_N = norm_width) %>% 
  dplyr::select(coordinates:width_without_N, contains("rpkm"), contains("rpm")) %>%  
  dplyr::arrange(-width_without_N) 
clusters_21dpp %>% 
  dplyr::select(coordinates, width_full, width_without_N, rpkm_21dpp = rpkm_WT.21dpp, rpm_21dpp = rpm_WT.21dpp) %T>% 
  readr::write_csv(., file.path(outpath, "piRNA_clusters.21dpp.rpkm_cutoff.1.20201231.csv"))


### join and plot together
# join tables
clusters_all <- 
  bind_rows(clusters_13dpp %>% dplyr::mutate(stage = "13dpp"), 
            clusters_21dpp %>% dplyr::mutate(stage = "21dpp")) %>% 
  dplyr::select(coordinates, stage, everything(.))

# create plot table
plot_tb <- 
  clusters_all %>% 
  dplyr::select(logfc_x = rpkm_log2FC.21dpp, logfc_y = rpkm_log2FC.13dpp, stage)

# set limits
axis_limits <- 13

# crosshair plot
cross_plot <-
  ggplot(plot_tb, aes(x = logfc_x, y = logfc_y, color = stage)) +
  geom_point(shape = 16, size = 3, alpha = 0.75) +
  scale_color_manual(values = c("13dpp" = "red3", "21dpp" = "blue")) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("log2FC_21dpp") +
  ylab("log2FC_13dpp") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, "rpkm_cutoff", rpkm_cutoff, "both_stages", "merged_clusters", "jpg", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)


