### INFO: 
### DATE: Mon May 04 01:07:48 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/Freimer_microarray")

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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# microarray results path
array_path <- file.path(inpath, "Freimer_2017.diff_exp.CrePos_vs_CreNeg.miR-15a.csv")

# ERCC estimate number of molecules
ercc_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Wu_2019_unpub_GSE133748/Analysis/ERCC_spike_expression",
                       "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.ERCC.Wu_2019.n_molecules.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read array results
array_tb <- readr::read_csv(array_path)

# read ERCC results
errc_tb <- readr::read_csv(ercc_path)

######################################################## MAIN CODE
### tidy data
# calculte mean number of molecules
ercc_tb_mean <- 
  errc_tb %>% 
  tidyr::pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "n_molecules") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(n_molecules = mean(n_molecules) %>% round(., 3))

# get miroarray expression values, join with ERCC molecule number estimates
array_tb_tidy <- 
  array_tb %>% 
  dplyr::filter(!is.na(gene_id)) %>% 
  dplyr::select(gene_id, gene_name, array_exprs = log2_mean_exp, padj, log2FoldChange) %>% 
  dplyr::mutate(array_exprs = 2^array_exprs) %>% 
  dplyr::left_join(., ercc_tb_mean, by = "gene_id") %>% 
  dplyr::mutate(array_exprs = log2(array_exprs),
                n_molecules = log2(n_molecules))

# build linear regression model
linearMod <- lm(array_exprs ~ n_molecules, data = array_tb_tidy)  # build linear regression model on full data
summary(linearMod)

# get annoation for significantly expressed genes
array_tb_annotation <- 
  array_tb_tidy %>% 
  dplyr::filter(padj <= 0.05,
                log2FoldChange < 0)

# plot 
scatter_plot <- 
  ggplot() + 
  geom_point(data = array_tb_tidy, aes(x = array_exprs, y = n_molecules), size = 1) +
  # coord_cartesian(xlim = c(0, 6),
  #                 ylim = c(0, 6)) +
  # geom_label_repel(data = array_tb_annotation, 
  #                  mapping = aes(x = array_exprs, y = n_molecules, label = gene_name)) +
  geom_point(data = array_tb_annotation, 
             mapping = aes(x = array_exprs, y = n_molecules), color = "red3", size = 3) +
  geom_abline(intercept = 0, color = "black", size = 1.5) +
  xlab("array mean expression (Freimer)") + 
  ylab("estimated no. molecules (ERCC)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "scatter_plot.array_expression.n_molecules.log2.png"), width = 10, height = 10)
