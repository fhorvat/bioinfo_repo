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


### 13dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.13dpp <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

# results path
results_path.13dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)


### 8.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.8.5dpp <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

# results path
results_path.8.5dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)

######################################################## READ DATA
### 8.5dpp 
# read fpkm table
fpkm_tb.8.5dpp <- readr::read_csv(fpkm_tb_path.8.5dpp)

# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp) %>% as_tibble(.)


### 13dpp 
# read fpkm table
fpkm_tb.13dpp <- readr::read_csv(fpkm_tb_path.13dpp)

# read results table 
results.13dpp <- openxlsx::read.xlsx(results_path.13dpp) %>% as_tibble(.)


######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# set FPKM cutoff
fpkm_cutoff <- 10

### plot crosshair plot - 8.5dpp vs. 13dpp
fpkm_tb_plot <- 
  left_join(fpkm_tb.8.5dpp, fpkm_tb.13dpp, by = "gene_id") %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > fpkm_cutoff)) %>% 
  dplyr::mutate(log2FC_8.5dpp = log2(KO_8.5dpp / WT_8.5dpp), 
                log2FC_13dpp = log2(KO_13dpp / WT_13dpp), 
                log2FC_21dpp = log2(KO_21dpp / WT_21dpp)) %>% 
  dplyr::select(gene_id, log2FC_8.5dpp, log2FC_13dpp, log2FC_21dpp) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  dplyr::left_join(., results.8.5dpp %>% dplyr::select(gene_id, log2FoldChange), by = "gene_id") %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
  dplyr::arrange(regulation)
  
# get limits
axis_limits <-
  c(fpkm_tb_plot$log2FC_8.5dpp, fpkm_tb_plot$log2FC_13dpp) %>%
  replace((is.infinite(.) | is.na(.)), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

# crosshair plot
cross_plot <-
  ggplot(fpkm_tb_plot, aes(x = log2FC_8.5dpp, y = log2FC_13dpp, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2 (Mov10l1 KO / WT) FPKM - 8.5dpp")) +
  ylab(str_c("log2 (Mov10l1 KO / WT) FPKM - 13dpp")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT.8.5dpp_vs_13dpp.log2_ratio.crosshair", "any_WT", str_c(fpkm_cutoff, "_rpkm_cutoff"), "jpg", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

