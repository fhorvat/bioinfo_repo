### INFO: draw Venn diagrams of significantly differentially expressed genes
### DATE: Sun Mar 11 05:42:56 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/results")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# read significantly up- and down-regulated genes
signif_genes <- list.files(outpath, pattern = "*.signif.csv")
signif_genes <- 
  lapply(signif_genes, readr::read_csv) %>% 
  magrittr::set_names(str_extract(signif_genes, "GV|MII|1C"))

######################################################## MAIN CODE
### create list of gene_id
# down-regulated
down_genes <- lapply(signif_genes, function(x){x %>% 
    dplyr::filter(log2FoldChange < 0) %$% gene_id})

# up-regulated
up_genes <- lapply(signif_genes, function(x){x %>% 
    dplyr::filter(log2FoldChange > 0) %$% gene_id})

### plot
# down-regulated
png(file = file.path(outpath, "Venn.CNOT6L.significant.downregulated.KO_vs_WT.GRCm38.89.png"), 
    width = 1000, height = 1000, units = "px", type = "cairo")
down_plot <- venn.diagram(x = down_genes, 
                          filename = NULL, 
                          fill = c("blue", "red", "green"), 
                          alpha = c(0.5, 0.5, 0.5), 
                          cex = 4,
                          cat.cex = 2,
                          category.names = names(down_genes), 
                          main = "Downregulated", 
                          main.cex = 2)
grid.draw(down_plot)
dev.off()

# up-regulated 
png(file = file.path(outpath, "Venn.CNOT6L.significant.upregulated.KO_vs_WT.GRCm38.89.png"), 
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = up_genes, 
                        filename = NULL, 
                        fill = c("blue", "red", "green"), 
                        alpha = c(0.5, 0.5, 0.5), 
                        cex = 4,
                        cat.cex = 2,
                        category.names = names(up_genes), 
                        main = "Upregulated", 
                        main.cex = 2)
grid.draw(up_plot)
dev.off()

