### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters/compare_clusters_with_RNAseq")

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

# small RNA-seq RPM table path
rpm_smallRNA_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Analysis/expression.1k_window/MesAur1.1k_windows.FPM_mean.csv"

# RNA-seq RPM table path
rpm_RNA_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression.1k_window/MesAur1.1k_windows.FPM_mean.csv"

# normalized width path
norm_width_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters/count_N_in_genome/MesAur1.1k_windows.norm_width.csv"

######################################################## READ DATA
# read normalized width table
norm_width_tb <- readr::read_csv(norm_width_path)

# read small RNA-seq RPM table
rpm_smallRNA_tb <- readr::read_csv(rpm_smallRNA_tb_path)

# read RNA-seq RPM table
rpm_RNA_tb <- readr::read_csv(rpm_RNA_tb_path)

######################################################## MAIN CODE
### clean tables
# small RNA-seq RPM to RPKM
rpm_smallRNA_tb %<>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::select(gene_id, coordinates, KO_13dpp_smallRNA = KO_13dpp, WT_13dpp_smallRNA = WT_13dpp) %>% 
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  tidyr::pivot_longer(cols = -c(gene_id, coordinates, width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(gene_id, coordinates), names_from = "sample_id", values_from = c("rpkm")) %>% 
  dplyr::mutate(log2FC_13dpp_smallRNA = log2(KO_13dpp_smallRNA / WT_13dpp_smallRNA))

# RNA-seq RPM to RPKM
rpm_RNA_tb %<>%  
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  dplyr::select(gene_id, coordinates, KO_13dpp_RNA = KO_13dpp, WT_13dpp_RNA = WT_13dpp) %>% 
  dplyr::left_join(., norm_width_tb, by = "gene_id") %>% 
  tidyr::pivot_longer(cols = -c(gene_id, coordinates, width), values_to = "rpm", names_to = "sample_id") %>% 
  dplyr::mutate(rpkm = round((rpm / (width / 1000)), 3)) %>% 
  tidyr::pivot_wider(id_cols = c(gene_id, coordinates), names_from = "sample_id", values_from = c("rpkm")) %>% 
  dplyr::mutate(log2FC_13dpp_RNA = log2(KO_13dpp_RNA / WT_13dpp_RNA))

### set parameters
# set name
table_name <- "MesAur1.1k_windows"

# set RPKM cutoff
rpkm_cutoff_list <- c(1, 5, 10)

### plot crosshair plot - 13dpp logFC smallRNA vs. 13dpp logFC RNA
purrr::map(rpkm_cutoff_list, function(rpkm_cutoff){
  
  rpkm_cutoff <- 5
  
  # get values
  rpm_tb_plot <- 
    left_join(rpm_smallRNA_tb, rpm_RNA_tb, by = c("gene_id", "coordinates")) %>% 
    dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > rpkm_cutoff)) %>% 
    dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0))) %>% 
    dplyr::select(coordinates, log2FC_13dpp_smallRNA, log2FC_13dpp_RNA) 
  
  # get limits
  axis_limits <-
    c(rpm_tb_plot$log2FC_13dpp_smallRNA, rpm_tb_plot$log2FC_13dpp_RNA) %>%
    replace(is.infinite(.), 0) %>% 
    abs(.) %>%
    max(.) %>%
    ceiling(.)
  
  # save limited clusters
  rpm_tb_plot %<>% 
    dplyr::filter(log2FC_13dpp_smallRNA > 0, log2FC_13dpp_RNA > 2) %>% 
    dplyr::filter_at(.vars = vars(contains("log2FC")), .vars_predicate = all_vars(!is.infinite(.))) %T>%
    readr::write_csv(., file.path(outpath, str_c(table_name, "Mov10l_KO_WT.13dpp.small_RNAseq_vs_RNAseq.log2_ratio.crosshair", "any_WT", 
                                                 str_c(rpkm_cutoff, "_rpkm_cutoff.limited"), "csv", sep = ".")))
  
  # crosshair plot
  cross_plot <-
    ggplot(rpm_tb_plot, aes(x = log2FC_13dpp_RNA, y = log2FC_13dpp_smallRNA)) +
    geom_point(shape = 16, size = 3, color = "gray30", alpha = 0.75) +
    scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
    scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(str_c("log2 (Mov10l1 KO / WT) RPKM in 13dpp - RNA-seq")) +
    ylab(str_c("log2 (Mov10l1 KO / WT) RPKM in 13dpp - small RNA-seq")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = 0.3),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT.13dpp.small_RNAseq_vs_RNAseq.log2_ratio.crosshair", "any_WT", 
                                             str_c(rpkm_cutoff, "_rpkm_cutoff.limited"), "jpg", sep = ".")),
         plot = cross_plot,
         width = 10, height = 10)
  
  # return
  return(rpkm_cutoff)
  
})
