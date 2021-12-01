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

### 8.5dpp piRNA expression
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")

# RPM table path
rpm_tb_path.8.5dpp.pirna <- list.files(expression_path, ".*\\.FPM_mean\\.csv$", full.names = T)

######################################################## READ DATA
### 8.5dpp 
# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp, sheet = "Mov10l1_KO_8.5dpp_vs_Mov10l1_WT") %>% as_tibble(.)

### 8.5dpp piRNA
# read fpkm table
rpm_tb.8.5dpp.pirna <- readr::read_csv(rpm_tb_path.8.5dpp.pirna)

######################################################## MAIN CODE
### clean tables
# get p-adjusted cut-off
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

# piRNA RPM
rpm_tb.8.5dpp.pirna %<>% 
  dplyr::select(gene_id, Mov10l1_WT.piRNA = Mov10l_WT_8.5, Mov10l1_KO.piRNA = Mov10l_KO_8.5)


### plot MA plot - 8.5dpp vs. 13dpp
# set name
table_name <- "MesAur1.RNA_seq_Mov10l_logFC_8.5dpp_vs_piRNA_8.5dpp_RPM"

# create table
ma_plot_tb <- 
  left_join(results.8.5dpp_tidy, rpm_tb.8.5dpp.pirna, by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(gene_id, log2FC_8.5dpp, Mov10l1_WT.piRNA, regulation) %>% 
  dplyr::arrange(regulation)

# get limits
y_axis_limit <-
  c(ma_plot_tb$log2FC_8.5dpp) %>%
  replace((is.infinite(.) | is.na(.)), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

x_axis_limit <- ma_plot_tb$Mov10l1_WT.piRNA %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.)
y_axis_limit <- 4

# crosshair plot
cross_plot <-
  ggplot(ma_plot_tb, aes(x = Mov10l1_WT.piRNA, y = log2FC_8.5dpp, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(trans = "log10",
                     limits = c(0.01, x_axis_limit),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-y_axis_limit, y_axis_limit), breaks = seq(-y_axis_limit, y_axis_limit, 2)) +
  xlab(str_c("piRNA Mov10l1 WT RPM - 8.5 dpp")) + 
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
                                           "MA_plot.DESeq2", 
                                           str_c("padj_", padj_cutoff), 
                                           "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

