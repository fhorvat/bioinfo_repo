### INFO: produce scatterplots of ratio of expression in Fugaku and CNOT6L data
### DATE: Sun Mar 11 00:36:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/fpkm")

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
library(GeneOverlap)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# set experiment name
experiment_name <- "Morgan2017"

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

# results path 
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# CNOT6L significantly diff. exp. genes
results_CNOT6L <- 
  list.files(file.path(results_path, "results"), "diffExp.CNOT6L.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# filter results table
results <- 
  results_CNOT6L %>% 
  dplyr::mutate(stage = str_extract(comparison, "1C|MII|GV")) %>% 
  dplyr::select(gene_id, log2FoldChange, padj, stage) %>% 
  dplyr::group_by(stage) %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::summarise(sign_count = n())

