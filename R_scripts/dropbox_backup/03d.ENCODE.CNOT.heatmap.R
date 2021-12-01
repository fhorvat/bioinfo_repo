### INFO: get expression in ENCODE mouse dataset
### DATE: Sun Jun 24 16:14:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/ENCODE_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(pheatmap)
library(RColorBrewer)

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

# summarizedOverlaps path
fpkm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Analysis/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.ENCODE_2014_mouse.FPKM.csv"

# CNOT genes info path
CNTO_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/CNOT.GRCm38.89.20180305.geneInfo.csv"

######################################################## READ DATA
# read FPKM 
fpkm_df <- readr::read_csv(file = fpkm_path)

# read info about CNOT genes
CNOT_df <- readr::read_csv(CNTO_info_path)

######################################################## MAIN CODE
# get FPKM of CNOT genes
fpkm_CNOT <- 
  fpkm_df %>% 
  dplyr::right_join(., CNOT_df %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::select(-gene_id) %>% 
  dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.)

# plot heatmap with annotation
pheatmap::pheatmap(mat = fpkm_CNOT,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   fontsize_row = 20, 
                   fontsize_col = 20,
                   col = colorRampPalette(brewer.pal(9, "Greys"))(20),
                   # breaks = seq(0, 100, by = 10),
                   cellwidth = 50, 
                   cellheight = 50, 
                   filename = file.path(outpath, "heatmap.CNOT_FPKM.ENCODE.GRCm38.89.png"),
                   height = 10,
                   width = 20)


### adjust CNOT8 in thymus - replace 30.0839118 with next smallest (17)
# get FPKM of CNOT genes
fpkm_CNOT <- 
  fpkm_df %>% 
  dplyr::right_join(., CNOT_df %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::select(-gene_id) %>% 
  dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>% 
  dplyr::mutate(thymus = replace(thymus, gene_name == "Cnot8", 17)) %>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.)

# plot heatmap with annotation
pheatmap::pheatmap(mat = fpkm_CNOT,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   fontsize_row = 20, 
                   fontsize_col = 20,
                   col = colorRampPalette(brewer.pal(9, "Greys"))(20),
                   # breaks = seq(0, 100, by = 10),
                   cellwidth = 50, 
                   cellheight = 50, 
                   filename = file.path(outpath, "heatmap.CNOT_FPKM.ENCODE.GRCm38.89.adjusted.png"),
                   height = 10,
                   width = 20)




