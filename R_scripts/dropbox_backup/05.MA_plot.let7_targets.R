### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/let7_targets")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(plotly)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### experiment
# set experiment name
experiment <- "Stein_2015_PLoSGenet_GSE57514"

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10_new")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis")

# results path
results_path <- file.path(analysis_path, "results")
  

### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.se.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.avgFPKM.csv$"), full.names = T)

# differential expression results path
diffExp_results_path <- list.files(path = results_path, pattern = str_c("diffExp.*", "ensembl.", ensembl_version, ".all.csv"), full.names = T)
  

### genome
# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# transcript info path
transcript_info_path <- list.files(genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.transcriptInfo.csv$"), full.names = T)


### targetScan
targetScan_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/let7_targets/TargetScan7.1__let-7-5p_miR-98-5p.predicted_targets.txt"

######################################################## READ DATA
# read differential expression path
diffExp_results <- readr::read_csv(diffExp_results_path)

# read transcript info
transcript_info <- readr::read_csv(transcript_info_path)

# read targetScan
targetScan_tb <- readr::read_delim(targetScan_path, delim = "\t")

######################################################## MAIN CODE
# get gene_id to targetScan table
targetScan_tb_tidy <- 
  targetScan_tb %>% 
  dplyr::select(gene_name.ts = `Target gene`, transcript_id = `Representative transcript`) %>% 
  dplyr::mutate(transcript_id = str_remove(transcript_id, "\\.[0-9]+$")) %>% 
  dplyr::left_join(., transcript_info, by = "transcript_id")
  
### DESeq2 output
# data for plot
plot_df <- 
  diffExp_results %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                significant = (padj < 0.1),
                let7_target = (gene_id %in% targetScan_tb_tidy$gene_id), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, !significant, "not_sign"),
                regulation = replace(regulation, let7_target, "let7_target"),
                regulation = factor(regulation, levels = c("not_sign", "up", "down", "let7_target"))) %>%
  dplyr::arrange(regulation)

# plot
ma_plot <- 
  ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  geom_point(size = 3, shape = 20) +
  scale_x_log10(limits = c(1e-01, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-5, 5),
                     breaks = c(-5:5)) +
  # scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "lightskyblue", let7_target = "black")) +
  scale_colour_manual(values = c(not_sign = "gray50", up = "gray50", down = "gray50", let7_target = "red3")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2FC") +
  ggtitle(str_c("Dicer", " KO vs. WT")) +
  theme_bw() +
  # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
  #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot", experiment, "KO_vs_WT", "ensembl", ensembl_version, "let7_targets.DESeq2.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)


### FPKMs
# data for plot
plot_df <- 
  diffExp_results %>%
  dplyr::select(Dicer_KO, Dicer_WT, padj, gene_id) %>% 
  dplyr::mutate(mean = ((Dicer_KO + Dicer_WT) / 2), 
                lfc = log2(Dicer_KO + 1) - log2(Dicer_WT + 1),
                padj = replace(padj, is.na(padj), 1), 
                significant = (padj < 0.1),
                let7_target = (gene_id %in% targetScan_tb_tidy$gene_id), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, !significant, "not_sign"),
                regulation = replace(regulation, let7_target, "let7_target"),
                regulation = factor(regulation, levels = c("not_sign", "up", "down", "let7_target"))) %>%
  dplyr::arrange(regulation)

# plot
ma_plot <- 
  ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  geom_point(size = 3, shape = 20) +
  scale_x_log10(limits = c(1e-01, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-5, 5),
                     breaks = c(-5:5)) +
  # scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "lightskyblue", let7_target = "black")) +
  scale_colour_manual(values = c(not_sign = "gray50", up = "gray50", down = "gray50", let7_target = "red3")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2FC") +
  ggtitle(str_c("Dicer", " KO vs. WT")) +
  theme_bw() +
  # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
  #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot", experiment, "KO_vs_WT", "ensembl", ensembl_version, "let7_targets.FPKMs.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)

