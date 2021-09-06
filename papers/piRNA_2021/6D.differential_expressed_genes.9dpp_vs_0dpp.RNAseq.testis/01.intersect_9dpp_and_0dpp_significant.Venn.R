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


### 0.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.0.5dpp.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded.all_samples")

# results path
results_path.0.5dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)
results_path.0.5dpp <- results_path.0.5dpp[!str_detect(results_path.0.5dpp, "old")]


### 8.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")

# results path
results_path.8.5dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)
results_path.8.5dpp <- results_path.8.5dpp[!str_detect(results_path.8.5dpp, "old")]


######################################################## READ DATA
# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp, sheet = "Mov10l1_KO_8.5dpp_vs_Mov10l1_WT") %>% as_tibble(.)
results.0.5dpp <- openxlsx::read.xlsx(results_path.0.5dpp, sheet = "Mov10l_KO_0.5dpp_vs_Mov10l_WT_0") %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

# create results list
results_list <- list(Mov10l1_KO_vs_WT_0dpp = results.0.5dpp, 
                     Mov10l1_KO_vs_WT_9dpp = results.8.5dpp)

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
png(filename = file.path(outpath, str_c(table_name, "Mov10l1_0dpp_vs_9dpp.Venn.upregulated", "png", sep = ".")), 
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
png(filename = file.path(outpath, str_c(table_name, "Mov10l1_0dpp_vs_9dpp.Venn.downregulated", "png", sep = ".")), 
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


