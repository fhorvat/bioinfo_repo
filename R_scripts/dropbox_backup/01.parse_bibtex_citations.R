### INFO: 
### DATE: Wed Oct 24 14:39:04 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/PhD_talks.2021")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(RefManageR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bibtex file path
bib_path <- file.path(inpath, "citations.bib")

######################################################## READ DATA
# read bibtex
bib <- RefManageR::ReadBib(bib_path)

######################################################## MAIN CODE
