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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### Mov10l1 KO path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3")

# results path
results_path.Mov10l1 <- list.files(expression_path, "protein_coding.*\\.all_results\\.xlsx$", full.names = T, recursive = T)
results_path.Mov10l1 <- results_path.Mov10l1[!str_detect(results_path.Mov10l1, "old")]


### Piwil1
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_MII_Piwil1.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")

# results path
results_path.Piwil1 <- list.files(expression_path, "protein_coding.*\\.all_results\\.xlsx$", full.names = T, recursive = T)
results_path.Piwil1 <- results_path.Piwil1[!str_detect(results_path.Piwil1, "old")]


######################################################## READ DATA
### Mov10l1 
# read results table 
results.Mov10l1 <- openxlsx::read.xlsx(results_path.Mov10l1, sheet = "Mov10l_KO_vs_Mov10l_WT") %>% as_tibble(.)

### Piwil1 
# read results table 
results.Piwil1 <- openxlsx::read.xlsx(results_path.Piwil1, sheet = "Piwil1_KO_vs_Piwil1_HET") %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# set padj value
padj_cutoff <- 0.01

# get significant genes
results.Piwil1.sign <- 
  results.Piwil1 %>%
  dplyr::filter(padj <= padj_cutoff)
# dplyr::filter_at(.vars = vars(matches("\\.FPKM")), any_vars(. > 1))

# clean tables
results.Piwil1_tidy <- 
  results.Piwil1 %>%
  dplyr::select(gene_id, log2FC_Piwil1 = log2FoldChange, padj, gene_name, gene_biotype) %>% 
  dplyr::mutate(regulation = ifelse(log2FC_Piwil1 > 0, "up", "down"), 
                regulation = replace(regulation, !(gene_id %in% results.Piwil1.sign$gene_id), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down")))

# clean tables
results.Mov10l1_tidy <- 
  results.Mov10l1 %>%
  dplyr::select(gene_id, log2FC_Mov10l1 = log2FoldChange)


### plot crosshair plot - Mov10l1 vs. Piwil1 logFC 
# join results to one table
crosshair_tb_plot <- 
  left_join(results.Piwil1_tidy, results.Mov10l1_tidy, by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::arrange(regulation)

# get limits
axis_limits <-
  c(crosshair_tb_plot$log2FC_Mov10l1, crosshair_tb_plot$log2FC_Piwil1) %>%
  replace((is.infinite(.) | is.na(.)), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

axis_limits <- 5

# crosshair plot
cross_plot <-
  ggplot(crosshair_tb_plot, aes(x = log2FC_Mov10l1, y = log2FC_Piwil1, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2FC DESeq2 Mov10l1 KO vs. Mov10l1 WT in GV")) +
  ylab(str_c("log2FC DESeq2 Piwil1 KO vs. Piwil1 HET in MII (Siomi)")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("Mov10l1_vs_Piwil1.log2_ratio.DESeq2.crosshair", 
                                           "Piwil1_sign",
                                           str_c("padj_", padj_cutoff), 
                                           "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)


### add gene names as labels
# filter data
plot_df_labels <-
  crosshair_tb_plot %>%
  dplyr::filter(regulation != "no") %>%
  dplyr::filter(log2FC_Mov10l1 > 1) %>%
  # dplyr::filter(abs(log2FC_Piwil1.KO_vs_HET) > 1) %>%
  # dplyr::filter(gene_id %in% mov10l1_results$gene_id) %>%
  dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

# add labels
cross_plot_labeled <-
  cross_plot +
  geom_text(data = plot_df_labels,
            aes(x = log2FC_Mov10l1, y = log2FC_Piwil1, label = gene_name),
            check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
            colour = "black", fontface = "plain")

# save plot
ggsave(filename = file.path(outpath, str_c("Mov10l1_vs_Piwil1.log2_ratio.DESeq2.crosshair", 
                                           "Piwil1_sign",
                                           str_c("padj_", padj_cutoff), 
                                           "labeled", 
                                           "png", sep = ".")),
       plot = cross_plot_labeled,
       width = 10, height = 10)



### linear model
# filter table
crosshair_tb_plot_filt <- 
  crosshair_tb_plot %>% 
  dplyr::filter_at(.vars = vars(starts_with("log2FC")), all_vars(!is.na(.) & !is.infinite(.)))

# build linear regression model - line has to go through origin so the zero term is added to the formula
linearMod <- lm(log2FC_Mov10l1 ~ 0 + log2FC_Piwil1, data = crosshair_tb_plot_filt)  # build linear regression model on full data
summary(linearMod)

# calculate predictions intervals
pred.int <- predict(linearMod, interval = "prediction")

# add to the table
crosshair_tb_plot_filt <- 
  cbind(crosshair_tb_plot_filt, pred.int) %>% 
  as_tibble(.)

# get annotation table
annotations <- tibble(xpos = Inf,
                      ypos = -Inf,
                      annotateText = str_c("R^2 adjusted ", round(summary(linearMod)$adj.r.squared, 4)))

# crosshair plot
cross_plot <-
  ggplot() +
  geom_point(crosshair_tb_plot_filt, 
             mapping = aes(x = log2FC_Mov10l1, y = log2FC_Piwil1, color = regulation, alpha = regulation),
             shape = 16, size = 3) +
  stat_smooth(data = crosshair_tb_plot_filt, 
              mapping = aes(x = log2FC_Mov10l1, y = log2FC_Piwil1), 
              method = lm) +
  geom_text(data = annotations, 
            mapping = aes(x = xpos, y = ypos, label = annotateText),
            colour = "black", fontface = "italic", size = 2.5,
            hjust = 1.03, vjust = -0.5) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2FC DESeq2 Mov10l1 KO vs. Mov10l1 WT in GV")) +
  ylab(str_c("log2FC DESeq2 Piwil1 KO vs. Piwil1 HET in MII (Siomi)")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("Mov10l1_vs_Piwil1.log2_ratio.DESeq2.crosshair", 
                                           "Piwil1_sign",
                                           "regression",
                                           str_c("padj_", padj_cutoff), 
                                           "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

