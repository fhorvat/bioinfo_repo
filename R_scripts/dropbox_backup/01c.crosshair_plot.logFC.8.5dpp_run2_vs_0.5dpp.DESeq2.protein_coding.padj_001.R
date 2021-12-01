### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### 8.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")

# results path
results_path.8.5dpp <- list.files(expression_path, "protein_coding.*\\.all_results\\.xlsx$", full.names = T, recursive = T)
results_path.8.5dpp <- results_path.8.5dpp[!str_detect(results_path.8.5dpp, "old")]


### 0.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.0.5dpp.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded.all_samples")

# results path
results_path.0.5dpp <- list.files(expression_path, "protein_coding.*\\.all_results\\.xlsx$", full.names = T, recursive = T)
results_path.0.5dpp <- results_path.0.5dpp[!str_detect(results_path.0.5dpp, "old")]


######################################################## READ DATA
### 8.5dpp 
# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp, sheet = "Mov10l1_KO_8.5dpp_vs_Mov10l1_WT") %>% as_tibble(.)

### 0.5dpp 
# read results table 
results.0.5dpp <- openxlsx::read.xlsx(results_path.0.5dpp, sheet = "Mov10l_KO_0.5dpp_vs_Mov10l_WT_0") %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# set padj value
padj_cutoff <- 0.01

# get significant genes
results.8.5dpp.sign <- 
  results.8.5dpp %>%
  dplyr::filter(padj <= padj_cutoff)
  # dplyr::filter_at(.vars = vars(matches("\\.FPKM")), any_vars(. > 1))

# clean tables
results.8.5dpp_tidy <- 
  results.8.5dpp %>%
  dplyr::select(gene_id, log2FC_8.5dpp = log2FoldChange, padj, gene_name, gene_biotype) %>% 
  dplyr::mutate(regulation = ifelse(log2FC_8.5dpp > 0, "up", "down"), 
                regulation = replace(regulation, !(gene_id %in% results.8.5dpp.sign$gene_id), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down")))

# clean tables
results.0.5dpp_tidy <- 
  results.0.5dpp %>%
  dplyr::select(gene_id, log2FC_0.5dpp = log2FoldChange)


### plot crosshair plot - 8.5dpp vs. 0.5dpp logFC 
# join results to one table
crosshair_tb_plot <- 
  left_join(results.8.5dpp_tidy, results.0.5dpp_tidy, by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::arrange(regulation)

# get limits
axis_limits <-
  c(crosshair_tb_plot$log2FC_8.5dpp, crosshair_tb_plot$log2FC_0.5dpp) %>%
  replace((is.infinite(.) | is.na(.)), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

axis_limits <- 4

# crosshair plot
cross_plot <-
  ggplot(crosshair_tb_plot, aes(x = log2FC_0.5dpp, y = log2FC_8.5dpp, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2FC DESeq2 Mov10l1 KO vs. Mov10l1 WT  - 0.5 dpp")) +
  ylab(str_c("log2FC DESeq2 Mov10l1 KO vs. Mov10l1 WT  - 8.5 dpp")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, 
                                           "Mov10l_KO_WT.8.5dpp_run2_vs_0.5dpp.log2_ratio.DESeq2.crosshair", 
                                           str_c("padj_", padj_cutoff), 
                                           "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

