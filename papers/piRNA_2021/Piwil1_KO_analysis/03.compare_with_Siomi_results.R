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
library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/RefSeq"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = ".*UCSCseqnames.geneInfo.csv$", full.names = T)

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# Piwil1 expression path
piwil1_path <- file.path(base_path, "hamster_MII_Piwil1.RNAseq/Analysis/expression.RefSeq")
piwil1_fpkm_path <- list.files(piwil1_path, ".*\\.FPKM_mean\\.csv$", full.names = T)
piwil1_results_path <- list.files(piwil1_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)

# Siomi's results path
siomis_results_path <- file.path(piwil1_path, "L1_MIIOo_DEG.xlsx")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read results
piwil1_results <- openxlsx::read.xlsx(piwil1_results_path) %>% as_tibble(.)
siomis_results <- openxlsx::read.xlsx(siomis_results_path) %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# get Siomi's significant results
siomis_results_sign <- 
  siomis_results %>% 
  dplyr::filter(estimatedDEG == 1) %>% 
  dplyr::select(gene_id, log2FoldChange = m.value)

# put to list
results_list <- list(DESeq2_results = piwil1_results, EdgeR_results = siomis_results_sign)


## get list of upregulated genes
upregulated_list <- invisible(purrr::map(names(results_list), function(result){
  
  # get significant results
  results_df_sign <- 
    results_list[[result]] %>% 
    dplyr::filter(log2FoldChange > 0) %$%
    gene_id
  
})) %>% 
  set_names(., names(results_list))

# plot 
png(filename = file.path(outpath, str_c(table_name, "Piwil1_KO_HET.DESeq2_vs_EdgeR.Venn.upregulated", "png", sep = ".")), 
    width = 2000, height = 2000)
venn.plot <- venn.diagram(upregulated_list, 
                          NULL, 
                          fill = gg_color_hue(length(upregulated_list)), 
                          alpha = rep(0.4, length(upregulated_list)), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 1.5,
                          category.names = str_replace(names(upregulated_list), ",", "\n vs. \n"),
                          main.cex = 2)
grid.newpage()
venn.plot
grid.draw(venn.plot)
dev.off()


## get list of downregulated genes
downregulated_list <- invisible(purrr::map(names(results_list), function(result){
  
  # get significant results
  results_df_sign <- 
    results_list[[result]] %>% 
    dplyr::filter(log2FoldChange < 0) %$%
    gene_id
  
})) %>% 
  set_names(., names(results_list))

# plot 
png(filename = file.path(outpath, str_c(table_name, "Piwil1_KO_HET.DESeq2_vs_EdgeR.Venn.downregulated", "png", sep = ".")), 
    width = 2000, height = 2000)
venn.plot <- venn.diagram(downregulated_list, 
                          NULL, 
                          fill = gg_color_hue(length(downregulated_list)), 
                          alpha = rep(0.4, length(downregulated_list)), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 1.5,
                          category.names = str_replace(names(downregulated_list), ",", "\n vs. \n"),
                          main.cex = 2)
grid.newpage()
venn.plot
grid.draw(venn.plot)
dev.off()

