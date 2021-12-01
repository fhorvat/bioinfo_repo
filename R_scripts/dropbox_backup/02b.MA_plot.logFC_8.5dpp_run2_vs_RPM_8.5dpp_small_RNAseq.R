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

# FPKM table path
fpkm_tb_path.8.5dpp <- list.files(expression_path, ".*\\.FPKM_long\\.csv$", full.names = T)

# results path
results_path.8.5dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)


### 8.5dpp piRNA expression
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")

# FPKM table path
fpkm_tb_path.8.5dpp.pirna <- list.files(expression_path, ".*\\.FPM_mean\\.csv$", full.names = T)

######################################################## READ DATA
### 8.5dpp RNA-seq
# read fpkm table
fpkm_tb.8.5dpp <- readr::read_csv(fpkm_tb_path.8.5dpp)

# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp) %>% as_tibble(.)


### 8.5dpp piRNA
# read fpkm table
fpkm_tb.8.5dpp.pirna <- readr::read_csv(fpkm_tb_path.8.5dpp.pirna)

######################################################## MAIN CODE
### clean tables
# RNA-seq FPKM
fpkm_tb.8.5dpp %<>% 
  dplyr::filter(sample_id != "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE") %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_KO|Mov10l1_WT"), 
                age = "8.5dpp") %>% 
  tidyr::unite(genotype_age, genotype, age) %>% 
  dplyr::group_by(gene_id, genotype_age) %>% 
  dplyr::summarise(fpkm = mean(fpkm)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = gene_id, names_from = genotype_age, values_from = fpkm) %>% 
  dplyr::select(gene_id, Mov10l1_WT = Mov10l1_WT_8.5dpp, Mov10l1_KO = Mov10l1_KO_8.5dpp)

# piRNA FPKM
fpkm_tb.8.5dpp.pirna %<>% 
  dplyr::select(gene_id, Mov10l1_WT.piRNA = Mov10l_WT_8.5, Mov10l1_KO.piRNA = Mov10l_KO_8.5, gene_name, gene_biotype)

# set name
table_name <- "MesAur1.RNA_seq_Mov10l_logFC_8.5dpp_vs_piRNA_8.5dpp_RPM"

# # set FPKM cutoff
# fpkm_cutoff <- 10

### plot crosshair plot - 8.5dpp vs. 13dpp
fpkm_tb_plot <- 
  left_join(fpkm_tb.8.5dpp, fpkm_tb.8.5dpp.pirna, by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>%
  # dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > fpkm_cutoff)) %>% 
  dplyr::mutate(log2FC_8.5dpp = log2(Mov10l1_KO / Mov10l1_WT)) %>% 
  dplyr::select(gene_id, log2FC_8.5dpp, Mov10l1_WT.piRNA) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  dplyr::left_join(., results.8.5dpp %>% dplyr::select(gene_id, log2FoldChange), by = "gene_id") %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
  dplyr::arrange(regulation)

# get limits
y_axis_limit <-
  c(fpkm_tb_plot$log2FC_8.5dpp) %>%
  replace((is.infinite(.) | is.na(.)), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

x_axis_limit <- fpkm_tb_plot$Mov10l1_WT.piRNA %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.)

# crosshair plot
cross_plot <-
  ggplot(fpkm_tb_plot, aes(x = Mov10l1_WT.piRNA, y = log2FC_8.5dpp, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(trans = "log10",
                     limits = c(0.01, x_axis_limit),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-y_axis_limit, y_axis_limit), breaks = seq(-y_axis_limit, y_axis_limit, 2)) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) +
  xlab(str_c("piRNA Mov10l1 WT RPM - 8.5 dpp")) + 
  ylab(str_c("log2 (Mov10l1 KO / WT) FPKM - 8.5 dpp")) +
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
                                           # "any_WT", str_c(fpkm_cutoff, "_rpkm_cutoff"), 
                                           "MA_plot", 
                                           "jpg", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

