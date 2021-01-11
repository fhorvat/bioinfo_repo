### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse")

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

library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# chosen genes path
genes_path <- list.files(inpath, "DevCell_genes_200620.FPKM.Gan.hamster_testis.PS.csv", full.names = T)

# mouse expression path
mouse_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Gan_2013_NatCommun_GSE35005/Analysis/expression"
mouse_tb_path <- list.files(mouse_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

# hamster expression path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression.added_RefSeq"
hamster_tb_path <- list.files(hamster_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

######################################################## READ DATA
# read genes table
genes_tb <- readr::read_csv(genes_path)

# read mouse expression
mouse_tb <- readr::read_csv(mouse_tb_path)

# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

######################################################## MAIN CODE
### mouse Gan spermatogenesis heatmap
# add mouse expression values
genes_mouse <- 
  genes_tb %>% 
  dplyr::select(gene_name, gene_id = mouse_gene_id) %>% 
  dplyr::left_join(., mouse_tb %>% dplyr::select(gene_id, primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC), by = "gene_id") %>% 
  # dplyr::filter(!(gene_name %in% c("Piwil3", "Tex101"))) %>% 
  # dplyr::filter(!is.na(gene_id)) %>% 
  dplyr::mutate(gene_name = factor(gene_name, levels = gene_name)) %>% 
  dplyr::select(-gene_id) %>%
  tidyr::pivot_longer(., cols = -gene_name, values_to = "FPKM", names_to = "stage") %>%
  dplyr::mutate(stage = factor(stage, levels = c("primitive_SG_A", "SG_A", "SG_B", "leptotene_SC", "pachytene_SC"))) %>% # max FPKM = 142.198, Sycp3
  dplyr::select(gene_name, stage, FPKM) %>% 
  dplyr::mutate(relative_FPKM = replace(FPKM, FPKM > 100, 100))

# write
genes_mouse %>% 
  dplyr::select(-relative_FPKM) %>% 
  tidyr::pivot_wider(id_cols = gene_name, values_from = "FPKM", names_from = "stage") %T>% 
  readr::write_csv(., file.path(outpath, "DevCell_genes_200620.Gan_FPKM.heatmap.csv"))

# plot using ggplot2
heat_plot <- 
  ggplot(genes_mouse, aes(x = stage, y = gene_name)) + 
  geom_tile(aes(fill = relative_FPKM), colour = "grey45") + 
  coord_equal() + 
  scale_fill_gradient(low = "white", high = "black", na.value = "white") +
  theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
        axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "bottom", 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x = "", 
       y = "", 
       fill = "") + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(genes_mouse$gene_name))) 

# save
ggsave(filename = file.path(outpath, str_c("DevCell_genes_200620.Gan_FPKM.heatmap.png", sep = ".")), plot = heat_plot, width = 10, height = 10)


### hamster testis Mov10l1 FC heatmap
# add hamster expression values
genes_hamster <- 
  genes_tb %>% 
  dplyr::select(gene_name, gene_id = hamster_gene_id) %>% 
  dplyr::left_join(., hamster_tb %>% dplyr::select(gene_id, Mov10l_WT_13dpp, Mov10l_KO_13dpp, Mov10l_WT_21dpp, Mov10l_KO_21dpp), by = "gene_id") %>% 
  # dplyr::filter(!(gene_name %in% c("Piwil3", "Tex101"))) %>% 
  # dplyr::filter(!is.na(gene_id)) %>% 
  dplyr::mutate(gene_name = factor(gene_name, levels = gene_name)) %>% 
  dplyr::mutate(ratio_KO_WT.13dpp = Mov10l_KO_13dpp / Mov10l_WT_13dpp, 
                ratio_KO_WT.21dpp = Mov10l_KO_21dpp / Mov10l_WT_21dpp) %>% 
  dplyr::select(gene_name, ratio_KO_WT.13dpp, ratio_KO_WT.21dpp) %>%
  tidyr::pivot_longer(., cols = -gene_name, values_to = "expression_ratio", names_to = "stage") %>%
  dplyr::mutate(stage = factor(stage, levels = c("ratio_KO_WT.13dpp", "ratio_KO_WT.21dpp"))) %>% # max FPKM = 371.182, Tex101
  dplyr::select(gene_name, stage, expression_ratio)

# write
genes_hamster %>% 
  tidyr::pivot_wider(id_cols = gene_name, values_from = "expression_ratio", names_from = "stage") %T>% 
  readr::write_csv(., file.path(outpath, "DevCell_genes_200620.Mov10l1_KO.relative_change.heatmap.csv"))

# plot using ggplot2
heat_plot <- 
  ggplot(genes_hamster, aes(x = stage, y = gene_name)) + 
  geom_tile(aes(fill = expression_ratio), colour = "grey45") + 
  coord_equal() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 1) +
  theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
        axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "bottom", 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x = "", 
       y = "", 
       fill = "") + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(genes_mouse$gene_name))) 

# save
ggsave(filename = file.path(outpath, str_c("DevCell_genes_200620.Mov10l1_KO.relative_change.heatmap.png", sep = ".")), plot = heat_plot, width = 10, height = 10)

