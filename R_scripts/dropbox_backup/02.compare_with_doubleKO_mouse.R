### INFO: 
### DATE: Fri Sep 27 22:22:13 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression")

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

# significant hamster Mov10l KO vs. WT results path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression/results/diffExp.DESeq2.genotype.all_biotype.significant_results.xlsx"

# significant mouse DBL vs. SOM results path
mouse_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Analysis/expression/results/diffExp.DESeq2.genotype.all_biotype.significant_results.xlsx"

######################################################## READ DATA
# read hamster significant genes
hamster_tb <- openxlsx::read.xlsx(hamster_path, sheet = "Mov10l_KO_vs_Mov10l_WT")

# read mouse significant genes
mouse_tb <- openxlsx::read.xlsx(mouse_path, sheet = "DBL_vs_SOM")

######################################################## MAIN CODE
# overlap gene names
intersect(hamster_tb$gene_name, mouse_tb$gene_name)


