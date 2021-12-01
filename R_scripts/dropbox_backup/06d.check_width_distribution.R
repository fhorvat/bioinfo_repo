### INFO: 
### DATE: Fri Dec 11 17:05:38 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq/distance_to_retrotransposons")

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
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gene info path
gene_info_path <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.geneInfo.csv")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# result path
result_path <- "Analysis/expression.added_PIWIL3.stranded/results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq"

# table name
table_name <- "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx"

# experiments list
experiment_names <- c("hamster_testis_Mov10l.0.5dpp.RNAseq", 
                      "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq", 
                      "hamster_testis_Mov10l.RNAseq", 
                      "hamster_testis_Mov10l.RNAseq", 
                      "hamster_testis_Mov10l.RNAseq")

# set sheet names
sheet_names <- 
  c("Mov10l_KO_0.5dpp_vs_Mov10l_WT_0", 
    "Mov10l1_KO_8.5dpp_vs_Mov10l1_WT", 
    "Mov10l_KO_13dpp_vs_Mov10l_WT_13", 
    "Mov10l_KO_21dpp_vs_Mov10l_WT_21", 
    "Mov10l_KO_adult_vs_Mov10l_WT_ad") %>% 
  set_names(experiment_names)

# results path
results_path_list <- 
  file.path(base_path, experiment_names, result_path, table_name) %>% 
  set_names(experiment_names)


### add oocyte to paths
# get new path
new_path <- file.path(base_path, "hamster_oocyte_Mov10l.RNAseq", 
                      "Analysis/expression.added_PIWIL3/results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq", 
                      "protein_coding.diffExp.DESeq2.genotype.significant_results.xlsx")

# experiment names
experiment_names <- c(experiment_names, "hamster_oocyte_Mov10l.RNAseq")

# sheet names
sheet_names <- c(sheet_names, "Mov10l_KO_vs_Mov10l_WT") %>% set_names(., experiment_names)

# results path
results_path_list <- c(results_path_list, new_path) %>% set_names(., experiment_names)


### add MII to paths
# get new path
new_path <- file.path(base_path, "hamster_MII_Piwil1.RNAseq", 
                      "Analysis/expression.added_PIWIL3.stranded/results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq", 
                      "protein_coding.diffExp.DESeq2.genotype.significant_results.xlsx")

# experiment names
experiment_names <- c(experiment_names, "hamster_MII_Piwil1.RNAseq")

# sheet names
sheet_names <- c(sheet_names, "Piwil1_KO_vs_Piwil1_HET") %>% set_names(., experiment_names)

# results path
results_path_list <- c(results_path_list, new_path) %>% set_names(., experiment_names)

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read differently expressed genes
results_tb_all <- purrr::map(sheet_names, function(sheet_name){
  
  # get experiment name
  exp_name <- names(sheet_names[sheet_names == sheet_name])
  
  # read table
  results_tb <- 
    openxlsx::read.xlsx(xlsxFile = results_path_list[exp_name], sheet = sheet_name) %>% 
    as_tibble(.) %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
    dplyr::mutate(experiment_name = exp_name, 
                  comparison = sheet_name, 
                  regulation = ifelse(log2FoldChange > 0, "up", "down"),
                  regulation = factor(regulation, levels = c("no", "up", "down")), 
                  start = as.numeric(start), end = as.numeric(end),
                  gene_length = end - start + 1) %>% 
    dplyr::select(experiment_name, comparison, gene_id, regulation, gene_length)
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# prepare data for plot
plot_tb <- 
  results_tb_all %>% 
  dplyr::mutate(age = str_extract(comparison, "0.5dpp|8.5dpp|13dpp|21dpp|adult|oocyte") %>% 
                  str_replace(., "0.5dpp", "0dpp") %>% 
                  str_replace(., "8.5dpp", "9dpp")) %>% 
  dplyr::mutate(age = replace(age, experiment_name == "hamster_oocyte_Mov10l.RNAseq", "GV"),
                age = replace(age, experiment_name == "hamster_MII_Piwil1.RNAseq", "MII")) %>% 
  dplyr::mutate(age = factor(age, levels = c("0dpp", "9dpp", "13dpp", "21dpp", "adult", "GV", "MII")))

# plot as boxplots
boxplot_plot <-
  ggplot() +
  stat_boxplot(data = plot_tb, aes(x = age, y = gene_length, fill = regulation), 
               geom = "errorbar", size = 1.5) +
  geom_boxplot(data = plot_tb, aes(x = age, y = gene_length, fill = regulation), 
               outlier.colour = "black", size = 1.5) +
  # geom_jitter(data = plot_tb, aes(x = exp_id, y = gene_length, color = regulation), 
  #             alpha = 1, width = 0.1, height = 0, show.legend = F, shape = 16, size = 3) +
  scale_y_continuous(limits = c(0, 100000)) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())

# save plot
ggsave(plot = boxplot_plot,
       filename = file.path(outpath, str_c("diff_exp_genes.gene_length", "boxplot.png", sep = ".")),
       width = 10,
       height = 10)
