### INFO: 
### DATE: Tue Dec 08 12:54:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi/coverage")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# coverage paths
het_bw_path <- file.path(inpath, "s_GV_Mov10l1_HET_So811_r1.PE_bismark_bt2_pe.deduplicated.raw.bw")
ko_bw_path <- file.path(inpath, "s_GV_Mov10l1_KO_So821_r1.PE_bismark_bt2_pe.deduplicated.raw.bw")

######################################################## READ DATA
# get 

######################################################## MAIN CODE

