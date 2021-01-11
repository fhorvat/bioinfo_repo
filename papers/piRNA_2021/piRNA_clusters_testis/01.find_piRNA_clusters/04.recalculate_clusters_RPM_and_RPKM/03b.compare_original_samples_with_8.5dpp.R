### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/recalculate_clusters_RPM_and_RPKM")

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

library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# input table path
clusters_tb_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters_final_full_reduced_200730-for recalculation.xlsx")

# original sequencing clusters path
clusters_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.csv")

# resequenced clusters path
clusters_reseq_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.8.5dpp.csv")

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_tb_path) %>% 
  as_tibble(.)

# read original sequencing clusters expression
original_tb <- readr::read_csv(clusters_path)

# read resequnced samples clusters expression
reseq_tb <- readr::read_csv(clusters_reseq_path)

######################################################## MAIN CODE
# join expression of clusters in original and resequnced samples 
joined_tb <- 
  original_tb %>% 
  dplyr::select(coordinates, 
                rpkm_WT_13dpp, rpkm_KO_13dpp, rpkm_log2FC_13dpp, 
                rpm_WT_13dpp, rpm_KO_13dpp, rpm_log2FC_13dpp) %>% 
  dplyr::left_join(., reseq_tb %>% dplyr::select(coordinates, 
                                                 rpkm_WT_8.5dpp, rpkm_KO_8.5dpp, rpkm_log2FC_8.5dpp, 
                                                 rpm_WT_8.5dpp, rpm_KO_8.5dpp,  rpm_log2FC_8.5dpp), 
                   by = "coordinates")

### scatterplots
genotype <- "KO"
counts <- "RPKM"

# filter table for plot
plot_tb <- 
  joined_tb %>% 
  dplyr::select(contains(genotype) & contains(counts)) %>% 
  magrittr::set_colnames(c("original", "reseq"))

# get limits
axis_limits <-
  c(plot_tb$original, plot_tb$reseq) %>%
  replace(is.infinite(.), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

# plot and save
scatter_plot <- 
  ggplot(data = plot_tb, aes(x = original, y = reseq)) + 
  geom_abline(intercept = 0, color = "red", size = 1.5) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(0, axis_limits), breaks = seq(0, axis_limits, 10)) +
  scale_y_continuous(limits = c(0, axis_limits), breaks = seq(0, axis_limits, 10)) +
  xlab("13dpp") + 
  ylab("8.5dpp") +
  ggtitle(str_c("13dpp vs. 8.5dpp", genotype, counts, sep = " ")) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, 
                              str_c("MesAur1.1k_pachytene_clusters.200730", "scatter_plot.compare_13dpp_with_8.5dpp", counts, genotype, "png", sep = ".")), 
         width = 12, height = 12)


### crosshair plots
# plot RPKM
plot_tb <- 
  joined_tb %>% 
  dplyr::select(original = rpkm_log2FC_13dpp, resequecing = rpkm_log2FC_8.5dpp)

# get limits
axis_limits <-
  c(plot_tb$original, plot_tb$resequecing) %>%
  replace(is.infinite(.), 0) %>% 
  abs(.) %>%
  max(.) %>%
  ceiling(.)

# crosshair plot
crosshair_plot <-
  ggplot(plot_tb, aes(x = original, y = resequecing)) +
  geom_point(shape = 16, size = 3, color = "gray30", alpha = 0.75) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  xlab(str_c("log2 (13dpp Mov10l1 KO / WT) RPKM")) +
  ylab(str_c("log2 (8.5dpp Mov10l1 KO / WT) RPKM")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MesAur1.1k_pachytene_clusters.200730", 
                                           "scatter_plot.compare_13dpp_with_8.5dpp", 
                                           "log2_ratio", "png", sep = ".")),
       plot = crosshair_plot,
       width = 10, height = 10)




