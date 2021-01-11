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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# Mov10l1 expression path 
mov10l1_path <- file.path(base_path, "hamster_oocyte_Mov10l.RNAseq/Analysis/expression/results.ensembl.99.MesAur1.0.20200415.UCSCseqnames")
mov10l1_path <- list.files(mov10l1_path, ".*\\.significant_results\\.xlsx$", full.names = T)

# Piwil1 expression path
piwil1_path <- file.path(base_path, "hamster_MII_Piwil1.RNAseq/Analysis/expression/results.ensembl.99.MesAur1.0.20200415.UCSCseqnames")
piwil1_path <- list.files(piwil1_path, ".*\\.significant_results\\.xlsx$", full.names = T)


######################################################## READ DATA
# read fpkm tables
mov10l1_fpkm <- openxlsx::read.xlsx(mov10l1_path, "Mov10l_KO_vs_Mov10l_WT") %>% as_tibble(.)
piwil1_fpkm <- openxlsx::read.xlsx(piwil1_path) %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# create results list
results_list <- list(Mov10l_KO_vs_Mov10l_WT = mov10l1_fpkm, 
                     Piwil1_KO_vs_Piwil1_HET = piwil1_fpkm)

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
png(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT_vs_Piwil1_KO_HET.Venn.upregulated", "png", sep = ".")), 
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
png(filename = file.path(outpath, str_c(table_name, "Mov10l_KO_WT_vs_Piwil1_KO_HET.Venn.downregulated", "png", sep = ".")), 
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


