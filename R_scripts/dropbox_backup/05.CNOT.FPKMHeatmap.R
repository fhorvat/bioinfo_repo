### INFO: heatmap of gene expression
### DATE: Sun Mar 11 15:06:03 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

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
wideScreen()

######################################################## FUNCTIONS
# ratio of rows/columns
plotHeatmap <- function(heatmap_matrix, heatmap_name){
  
  # plot heatmap with annotation
  pheatmap::pheatmap(heatmap_matrix,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize_row = 20, 
                     fontsize_col = 20,
                     col = colorRampPalette(brewer.pal(9, "Greys"))(20),
                     # breaks = seq(0, 100, by = 10),
                     cellwidth = 50, 
                     cellheight = 50, 
                     filename = file.path(outpath, "results", str_c("heatmap.", heatmap_name, ".png")),
                     height = 10,
                     width = 10)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# CNOT in Fugaku's data
CNOT_Fugaku <- readr::read_csv(file.path(inpath, "CNOT.GRCm38.89.Fugaku.FPKM.csv"))

# CNOT in Graf_2014 cow data
CNOT_Graf <- readr::read_csv(file.path(inpath, "CNOT.UMD3.1.91.bosTau8.Graf_2014.FPKM.csv"))

# CNOT in Xue_2013 human data
CNOT_Xue <- readr::read_csv(file.path(inpath, "CNOT.GRCh38.91.20180312.UCSCnames.clean.Xue2013.avgFPKM.csv"))

# CNOT in Dang_2016 human data
CNOT_Dang <- readr::read_csv(file.path(inpath, "CNOT.GRCh38.91.20180312.UCSCnames.clean.Dang2016.avgFPKM.csv"))

# CNOT in Hendrickson_2017 human data
CNOT_Hendrickson <- readr::read_csv(file.path(inpath, "CNOT.GRCh38.91.20180312.UCSCnames.clean.Hendrickson2017.avgFPKM.csv"))

######################################################## MAIN CODE
# Fugaku
CNOT_Fugaku %>% 
  dplyr::select(matches(".WE$|gene_name")) %>% 
  magrittr::set_colnames(., str_replace(colnames(.), ".WE", "")) %>% 
  dplyr::select(gene_name, GV, MII, `1C`, `2C`, `4C`, Morula = Molura, Blast) %>%
  dplyr::filter(str_detect(gene_name, "Cnot")) %>%
  dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>%
  # dplyr::right_join(tibble(gene_name = c(str_c("Cnot", c(1:4, 6, "6l", 7:11)), "Pan2", "Pan3", "Parn")), by = "gene_name") %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %>% 
  plotHeatmap(heatmap_matrix = ., 
              heatmap_name = "CNOT_FPKM.Fugaku.GRCm38.89.extended")

# Graf
CNOT_Graf %>% 
  dplyr::select(gene_name, GV, MII, `4C`, `8C`, `16C`, Blast = blast) %>%
  # dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  # dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>% 
  dplyr::right_join(tibble(gene_name = c(str_c("Cnot", c(1:4, 6, "6l", 7:11)), "Pan2", "Pan3", "Parn")), by = "gene_name") %>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %>% 
  plotHeatmap(heatmap_matrix = ., 
              heatmap_name = "CNOT_FPKM.Graf.UMD3.1.91.extended")

# Xue
CNOT_Xue %>% 
  dplyr::select(gene_name, MII, `1C`, `2C_blastomere`, `4C_blastomere`, `8C_blastomere`, morula) %>%
  # dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  dplyr::mutate(`1C` = replace(`1C`, gene_name == "Cnot11", 105)) %>% 
  # dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>%
  dplyr::right_join(tibble(gene_name = c(str_c("Cnot", c(1:4, 6, "6l", 7:11)), "Pan2", "Pan3", "Parn")), by = "gene_name") %>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %T>% 
  plotHeatmap(heatmap_matrix = ., 
              heatmap_name = "CNOT_FPKM.Xue.GRCh38.91.limit_1C_CNOT11.extended")

# Dang
CNOT_Dang %>% 
  dplyr::select(gene_name, MII, `1C`, `2C`, `4C`, `8C`, morula, early_blast, middle_blast, late_blast, ICM, troph, ESC) %>%
  # dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  dplyr::mutate(`1C` = replace(`1C`, gene_name == "Cnot11", 63)) %>% 
  # dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>% 
  dplyr::right_join(tibble(gene_name = c(str_c("Cnot", c(1:4, 6, "6l", 7:11)), "Pan2", "Pan3", "Parn")), by = "gene_name") %>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %T>% 
  plotHeatmap(heatmap_matrix = ., 
              heatmap_name = "CNOT_FPKM.Dang.GRCh38.91.limit_1C_CNOT11.extended")

# Hendrickson
CNOT_Hendrickson %>% 
  dplyr::select(gene_name, GV, MI, MII, `1C`, `2_4_8C`, Morula, ICM, troph) %>%
  dplyr::filter(str_detect(gene_name, "Cnot")) %>% 
  # dplyr::mutate(`1C` = replace(`1C`, gene_name == "Cnot11", 63)) %>% 
  dplyr::right_join(tibble(gene_name = str_c("Cnot", c(1:4, 6, "6l", 7:11))), by = "gene_name") %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %T>% 
  plotHeatmap(heatmap_matrix = ., 
              heatmap_name = "CNOT_FPKM.Hendrickson.GRCh38.91.extended")
