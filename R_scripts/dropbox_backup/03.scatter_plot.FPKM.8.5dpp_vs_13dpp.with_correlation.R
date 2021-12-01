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
  # dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > fpkm_cutoff)) %>% 
  dplyr::mutate(log2FC_8.5dpp = log2(KO_8.5dpp / WT_8.5dpp), 
                log2FC_13dpp = log2(KO_13dpp / WT_13dpp), 
                log2FC_21dpp = log2(KO_21dpp / WT_21dpp)) %>% 
  dplyr::select(gene_id, WT_8.5dpp, WT_13dpp, log2FC_8.5dpp, log2FC_13dpp, log2FC_21dpp) %>%
  # dplyr::mutate(WT_8.5dpp = log2(WT_8.5dpp + 0.001), 
  #               WT_13dpp = log2(WT_13dpp + 0.001)) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
  dplyr::left_join(., results.8.5dpp %>% dplyr::select(gene_id, log2FoldChange), by = "gene_id") %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down")))
  # dplyr::arrange(regulation)


### linear model
# build linear regression model - line has to go through origin so the zero term is added to the formula
linearMod <- lm(WT_8.5dpp ~ 0 + WT_13dpp, data = fpkm_tb_plot)  # build linear regression model on full data
summary(linearMod)

# calculate predictions intervals
pred.int <- predict(linearMod, interval = "prediction")

# add to the table
fpkm_tb_plot <- 
  cbind(fpkm_tb_plot, pred.int) %>% 
  as_tibble(.)

# get annotation table
annotations <- tibble(xpos = Inf,
                      ypos = -Inf,
                      annotateText = str_c("R^2 adjusted ", round(summary(linearMod)$adj.r.squared, 4)))

# plot regression line + confidence intervals
ln_plot <- 
  ggplot() +
  geom_point(data = fpkm_tb_plot, aes(x = WT_8.5dpp, y = WT_13dpp)) +
  stat_smooth(data = fpkm_tb_plot, aes(x = WT_8.5dpp, y = WT_13dpp), method = lm) +
  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
            colour = "black", fontface = "italic", size = 2.5,
            hjust = 1.03, vjust = -0.5) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  xlab(str_c("Mov10l1 WT 8.5dpp FPKM")) +
  ylab(str_c("Mov10l1 WT 13dpp FPKM")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_WT.8.5dpp_vs_13dpp.FPKM.scatterplot", 
                                             "lim_FPKM_100", 
                                             "jpg", sep = ".")),
         height = 10, width = 10)

# scatter plot
scatter_plot <-
  ggplot(fpkm_tb_plot, aes(x = WT_8.5dpp, y = WT_13dpp)) +
  geom_point(shape = 16, size = 3) +
  stat_smooth(data = fpkm_tb_plot, aes(x = WT_8.5dpp, y = WT_13dpp), method = lm) +
  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
            colour = "black", fontface = "italic", size = 2.5,
            hjust = 1.03, vjust = -0.5) +
  # scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
  #                     values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  # scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab(str_c("Mov10l1 WT 8.5dpp FPKM")) +
  ylab(str_c("Mov10l1 WT 13dpp FPKM")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_WT.8.5dpp_vs_13dpp.FPKM.scatterplot.log", "jpg", sep = ".")),
       plot = scatter_plot,
       width = 10, height = 10)


