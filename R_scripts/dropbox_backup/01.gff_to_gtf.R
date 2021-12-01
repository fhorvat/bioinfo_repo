### INFO: 
### DATE: Wed Jul 14 19:14:27 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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

library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# gff path
gff_path <- file.path(inpath, "AriVul.gff")

######################################################## READ DATA
# read  gff
gff_gr <- rtracklayer::import(gff_path)

######################################################## MAIN CODE
# save
rtracklayer::export(gff_gr, file.path(outpath, "AriVul.gtf"))
