### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Green_2018_DevCell_GSE112393/Analysis/expression/reproduce_paper")

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

library(Seurat)
library(patchwork)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(getwd(), "..", "Seurat_analysis")

# raw counts tabel
count_tb_path <- file.path(inpath, "GSE112393_MergedAdultMouseST25_DGE.filtered.matrix.RDS")

######################################################## READ DATA
# read filtered count data
count_tb <- readRDS(count_tb_path)

######################################################## MAIN CODE
# change order of dataset in orig.ident
dataset <- 
  colnames(count_tb) %>% 
  str_remove(., "_.*") %>% 
  unique(.) %>% 
  .[c(1:16, 25, 17:24)]

### Seurat analysis
# create Seurat object
count_seurat <- CreateSeuratObject(counts = count_tb, project = "mouse_testis_scRNAseq")

# change order of the dataset in Seurat object
count_seurat@meta.data$orig.ident <- factor(count_seurat@meta.data$orig.ident, levels = dataset)
count_seurat@active.ident <- factor(count_seurat@active.ident, levels = dataset)

# add percent.mito to the datainfo
count_seurat[["percent.mito"]] <- PercentageFeatureSet(count_seurat, pattern = "^mt-")


# normalize data
count_norm <- NormalizeData(count_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


