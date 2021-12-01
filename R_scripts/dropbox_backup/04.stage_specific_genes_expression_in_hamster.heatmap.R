### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Gan_2013_NatCommun_GSE35005")

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
clusters_tb_path <- file.path(inpath, "Gan_2013_NatCommun_GSE35005.stage_specific_genes.hamster_gene_id.csv")

######################################################## READ DATA
# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

# read clusters table
clusters_tb <- readr::read_csv(clusters_tb_path)

######################################################## MAIN CODE
# filter clusters table
cluster_fpkm <- 
  clusters_tb %>% 
  dplyr::left_join(., hamster_tb, by = c("hamster_gene_id" = "gene_id")) %>% 
  dplyr::select(gene_id = hamster_gene_id, gene_name, cluster = comparison, mouse_logFC, Mov10l_WT_13dpp:Mov10l_WT_adult, Mov10l_KO_13dpp:Mov10l_WT_adult) %>% 
  dplyr::mutate(gene_id = make.unique(gene_id)) %>% 
  dplyr::filter(!is.na(cluster)) %>% 
  magrittr::set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::mutate(cluster = factor(cluster, levels = rev(c("primitive_SG_A", "SG_A", "SG_B", 
                                                     "leptotene_SC", "pachytene_SC", 
                                                     "round_ST", "elongative_ST")))) %>% 
  dplyr::arrange(cluster, desc(mouse_logFC))

### FPKM heatmap
# create matrix
cluster_matrix <-
  cluster_fpkm %>%
  dplyr::select(gene_id, WT_13dpp:KO_adult) %>%
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "gene_id") %>%
  as.matrix(.)

# log transform
cluster_matrix <- log2(cluster_matrix + 0.1)

# annotation data.frame
annotation_df <-
  cluster_fpkm %>%
  dplyr::select(gene_id, stage = cluster) %>%
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "gene_id")

# plot
pheatmap::pheatmap(cluster_matrix,
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   col = colorRampPalette(c("white", "black"))(101),
                   border_color =  NA,
                   annotation_row = annotation_df,
                   file = file.path(outpath, str_c("ensembl", ensembl_version, "Gan_2013.stage_specific_genes.hamster_testis.FPKM", "png", sep = ".")),
                   height = 15,
                   width = 10)
