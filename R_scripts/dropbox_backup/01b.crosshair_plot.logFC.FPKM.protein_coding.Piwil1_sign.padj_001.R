### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Piwil1_KO_analysis")

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

library(openxlsx)
library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# Mov10l1 expression path 
mov10l1_path <- file.path(base_path, "hamster_oocyte_Mov10l.RNAseq/Analysis/expression.added_PIWIL3")
mov10l1_fpkm_path <- list.files(mov10l1_path, ".*\\.FPKM_mean\\.csv$", full.names = T)
mov10l1_results_path <- list.files(mov10l1_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)
mov10l1_results_path <- mov10l1_results_path[!(str_detect(mov10l1_results_path, "old"))]

# Piwil1 expression path
piwil1_path <- file.path(base_path, "hamster_MII_Piwil1.RNAseq/Analysis/expression.added_PIWIL3.stranded")
piwil1_fpkm_path <- list.files(piwil1_path, ".*\\.FPKM_mean\\.csv$", full.names = T)
piwil1_results_path <- list.files(piwil1_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)
piwil1_results_path <- piwil1_results_path[!(str_detect(piwil1_results_path, "old"))]

######################################################## READ DATA
# read fpkm tables
mov10l1_fpkm <- readr::read_csv(mov10l1_fpkm_path)
piwil1_fpkm <- readr::read_csv(piwil1_fpkm_path)

# read results
mov10l1_results <- openxlsx::read.xlsx(mov10l1_results_path, "Mov10l_KO_vs_Mov10l_WT") %>% as_tibble(.)
piwil1_results <- openxlsx::read.xlsx(piwil1_results_path) %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# set FPKM cutoff
fpkm_cutoff <- 0

# set p-adjusted cutoff 
padj_cutoff <- 0.01

### plot crosshair plot - Mov10l1 KO vs. Piwil1 KO
# get logFC values in Mov10l1
mov10l1_lfc <- 
  mov10l1_fpkm %>% 
  dplyr::select(gene_id, gene_name, Mov10l_KO, Mov10l_WT, Mov10l_HET) %>% 
  dplyr::mutate(log2FC_Mov10l.KO_vs_WT = log2(Mov10l_KO / Mov10l_WT), 
                log2FC_Mov10l.KO_vs_HET = log2(Mov10l_KO / Mov10l_HET))

# get logFC values in Piwil1
piwil1_lfc <- 
  piwil1_fpkm %>% 
  dplyr::select(gene_id, Piwil1_KO, Piwil1_HET_MII = Piwil1_HET) %>% 
  dplyr::mutate(log2FC_Piwil1.KO_vs_HET = log2(Piwil1_KO / Piwil1_HET_MII))

# join 
fpkm_tb_plot <- 
  dplyr::left_join(mov10l1_lfc, piwil1_lfc, by = "gene_id") %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  dplyr::left_join(., piwil1_results %>% dplyr::select(gene_id, log2FoldChange), by = "gene_id") %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
  dplyr::arrange(regulation)

# get limits
axis_limits <-
  c(fpkm_tb_plot$log2FC_Mov10l.KO_vs_WT, fpkm_tb_plot$log2FC_Piwil1.KO_vs_HET) %>%
  replace(is.infinite(.), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

# crosshair plot
cross_plot <-
  ggplot(fpkm_tb_plot, aes(x = log2FC_Mov10l.KO_vs_WT, y = log2FC_Piwil1.KO_vs_HET, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2 (Mov10l1 KO / WT) FPKM")) +
  ylab(str_c("log2 (Piwil1 KO / HET) FPKM")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("Mov10l1_vs_Piwil1.log2_ratio.FPKM.crosshair", 
                                           "Piwil1_sign",
                                           str_c("padj_", padj_cutoff), 
                                           "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)



# ### add gene names as labels
# # filter data
# plot_df_labels <-
#   fpkm_tb_plot %>%
#   # dplyr::filter(regulation != "no") %>%
#   # dplyr::filter(abs(log2FC_Mov10l.KO_vs_WT) > 1) %>% 
#   # dplyr::filter(abs(log2FC_Piwil1.KO_vs_HET) > 1) %>% 
#   dplyr::filter(gene_id %in% mov10l1_results$gene_id) %>% 
#   dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
# 
# # add labels
# cross_plot_labeled <-
#   cross_plot +
#   geom_text_repel(data = plot_df_labels,
#                   aes(x = log2FC_Mov10l.KO_vs_WT, y = log2FC_Piwil1.KO_vs_HET, label = gene_name),
#                   size = 2, hjust = 0, vjust = 1.5,
#                   colour = "black", alpha = 1, fontface = "plain")
# 
# # save plot
# ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT_vs_Piwil1_KO_HET.log2_ratio.crosshair.labeled", "jpg", sep = ".")),
#        plot = cross_plot_labeled,
#        width = 10, height = 10)

