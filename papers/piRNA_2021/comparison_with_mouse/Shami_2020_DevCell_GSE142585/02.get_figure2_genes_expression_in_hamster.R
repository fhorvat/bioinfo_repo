### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Shami_2020_DevCell_GSE142585")

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

library(pheatmap)
library(RColorBrewer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# hamster expression path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression"
hamster_tb_path <- list.files(hamster_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

# table of consensus SPG clusters
clusters_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Shami_2020_DevCell_GSE142585"
clusters_tb_path <- file.path(clusters_path, "Shami_2020_DevCell_GSE142585.figure2_clusters.hamster_gene_id.csv")

######################################################## READ DATA
# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

# read clusters table
clusters_tb <- readr::read_csv(clusters_tb_path)

######################################################## MAIN CODE
# set FPKM cut-off filter
cutoff <- 10

# filter clusters table
cluster_fpkm <- 
  clusters_tb %>% 
  dplyr::left_join(., hamster_tb, by = c("hamster_gene_id" = "gene_id")) %>% 
  dplyr::select(gene_id = hamster_gene_id, gene_name, cluster = Cluster, Mov10l_WT_13dpp:Mov10l_WT_adult, Mov10l_KO_13dpp:Mov10l_WT_adult) %>% 
  dplyr::mutate(gene_id = make.unique(gene_id)) %>% 
  dplyr::filter(!is.na(cluster)) %>% 
  magrittr::set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::arrange(cluster) %>% 
  dplyr::mutate(cluster = str_c("cl_", cluster)) %>% 
  dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > cutoff))


### log2FC heatmap
# calculate log2FC
cluster_fc <-
  cluster_fpkm %>%
  dplyr::mutate(log2FC_13dpp = log2(KO_13dpp / WT_13dpp),
                log2FC_21dpp = log2(KO_21dpp / WT_21dpp),
                log2FC_adult = log2(KO_adult / WT_adult)) %>%
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>%
  dplyr::filter_at(.vars = vars(contains("log2FC")), .vars_predicate = all_vars(!is.infinite(.))) %>% 
  dplyr::arrange(cluster, desc(log2FC_21dpp))

# save
readr::write_csv(cluster_fc, file.path(outpath, str_c("ensembl", ensembl_version, "Shami_2020.figure2_clusters.hamster_testis.expression.csv", sep = ".")))

# create matrix
cluster_matrix <-
  cluster_fc %>%
  dplyr::select(gene_id, log2FC_13dpp:log2FC_adult) %>%
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "gene_id") %>%
  as.matrix(.)

# annotation data.frame
annotation_df <-
  cluster_fpkm %>%
  dplyr::select(gene_id, cluster) %>%
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "gene_id")

# plot
pheatmap::pheatmap(cluster_matrix,
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   col = colorRampPalette(c("red", "white", "green"))(3),
                   breaks = c(-16, -1, 1, 5),
                   annotation_row = annotation_df,
                   file = file.path(outpath, str_c("ensembl", ensembl_version, "Shami_2020.figure2_clusters.hamster_testis.log2FC.png", sep = ".")),
                   height = 15,
                   width = 10)


# ### FPKM heatmap
# # create matrix
# cluster_matrix <- 
#   cluster_fpkm %>% 
#   dplyr::select(gene_id, WT_13dpp:KO_adult) %>% 
#   as.data.frame(.) %>% 
#   tibble::column_to_rownames(., var = "gene_id") %>% 
#   as.matrix(.)
# 
# # cluster_matrix[cluster_matrix > 100] <- 100
# 
# # # log transform
# # cluster_matrix <- log2(cluster_matrix + 0.1)
# 
# # annotation data.frame
# annotation_df <-
#   cluster_fpkm %>%
#   dplyr::select(gene_id, cluster) %>% 
#   as.data.frame(.) %>% 
#   tibble::column_to_rownames(., var = "gene_id")
# 
# # plot
# pheatmap::pheatmap(cluster_matrix,
#                    cluster_rows = F,
#                    cluster_cols = F,
#                    show_rownames = F, 
#                    show_colnames = T, 
#                    col = colorRampPalette(c("white", "black"))(101),
#                    breaks = c(1:100, 10000), 
#                    annotation_row = annotation_df, 
#                    file = file.path(outpath, str_c("ensembl", ensembl_version, "Shami_2020.figure2_clusters.hamster_testis.FPKM", str_c(cutoff, "_cutoff"), ".png", sep = ".")),
#                    height = 15,
#                    width = 10)
# 
# 
