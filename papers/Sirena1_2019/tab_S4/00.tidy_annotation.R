### INFO: 
### DATE: Wed Nov 27 18:27:46 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Sirena1_pseudogenes_expression")

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

# coordinates table path
coords_path <- file.path(inpath, "2019 Ganesh Table S4 Elob pseudogenes.xlsx")

######################################################## READ DATA
# read coordinates
coords_tb <- openxlsx::read.xlsx(coords_path) %>% as_tibble(.)
  
######################################################## MAIN CODE
# get coordinates into SAF format required by featureCounts
coords_saf <- 
  coords_tb %>% 
  dplyr::select(GeneID = ID, coordinates) %>% 
  tidyr::separate(coordinates, into = c("Chr", "coordinates"), sep = ":") %>% 
  tidyr::separate(coordinates, into = c("Start", "End"), sep = "-") %>% 
  dplyr::mutate(Strand = "+") %>% 
  readr::write_delim(., path = file.path(outpath, "2019_Ganesh_Table_S4_Elob_pseudogenes.saf"), delim = "\t")



