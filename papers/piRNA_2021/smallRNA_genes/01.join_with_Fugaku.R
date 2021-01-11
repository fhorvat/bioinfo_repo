### INFO: 
### DATE: Fri Apr 24 10:20:01 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/smallRNA_genes")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# original table path
original_tb_path <- file.path(inpath, "smallRNA.perfect.21to23.all.mouse.200424 working Table S1 draft.xlsx")

# Fugaku FPKM path
fugaku_fpkm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.Fugaku.FPKM.csv"

######################################################## READ DATA
# read original table
original_tb <- openxlsx::read.xlsx(original_tb_path)

# read Fugaku's FPKM table
fpkm_tb <- readr::read_csv(fugaku_fpkm_path)

######################################################## MAIN CODE
# filter
fpkm_tb_GV <- 
  fpkm_tb %>% 
  dplyr::select(gene.id = gene_id, s_GV.WE.PE) %>% 
  dplyr::mutate(s_GV.WE.PE = round(s_GV.WE.PE, 3))

# add to original table
original_tb_added <- 
  original_tb %>% 
  as_tibble(.) %>% 
  left_join(., fpkm_tb_GV, by = "gene.id")

# write table
readr::write_csv(original_tb_added, file.path(outpath, "smallRNA.perfect.21to23.all.mouse.200424 working Table S1 draft.added_Fugaku.csv"))