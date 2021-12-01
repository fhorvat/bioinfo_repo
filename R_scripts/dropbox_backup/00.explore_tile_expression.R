### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/RNAseq")

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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.1k_window")

# RPM table path
rpm_tb_path <- list.files(expression_path, ".*\\.FPM_mean\\.csv$", full.names = T)

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# sample table path
sample_tb_path <- list.files(documentation_path, ".*\\.sampleTable\\.csv$", full.names = T)

######################################################## READ DATA
# read RPM table
rpm_tb <- readr::read_csv(rpm_tb_path)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.1k_windows"

### plot crosshair plot - 13dpp vs. 21dpp
# get values
rpm_tb_plot <- 
  rpm_tb %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > 10)) %>% 
  dplyr::mutate(log2FC_13dpp = log2(KO_13dpp / WT_13dpp), 
                log2FC_21dpp = log2(KO_21dpp / WT_21dpp)) %>% 
  dplyr::select(coordinates, log2FC_13dpp, log2FC_21dpp) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0)))

# get limits
axis_limits <-
  c(rpm_tb_plot$log2FC_13dpp, rpm_tb_plot$log2FC_21dpp) %>%
  replace(is.infinite(.), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

# crosshair plot
cross_plot <-
  ggplot(rpm_tb_plot, aes(x = log2FC_21dpp, y = log2FC_13dpp)) +
  geom_point(shape = 16, size = 3, color = "gray30", alpha = 0.75) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2 (Mov10l1 KO / WT) RPM - 21dpp")) +
  ylab(str_c("log2 (Mov10l1 KO / WT) RPM - 13dpp")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT.13dpp_vs_21dpp.log2_ratio.crosshair", "any_WT", "10_rpkm_cutoff", "jpg", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

