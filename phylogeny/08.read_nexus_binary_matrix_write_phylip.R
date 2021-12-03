### INFO: 
### DATE: Tue May 05 19:14:57 2020
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

library(ape)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# Supertree binary matrix in nexus format path
binmat_path <- list.files(inpath, "\\.nex$", full.names = T)

######################################################## READ DATA
# read Supertree binary matrix in nexus format
binmat_nexus <- ape::read.nexus.data(binmat_path)

######################################################## MAIN CODE
# get name
binmat_name <- 
  binmat_path %>% 
  basename(.) %>% 
  str_remove(., "\\.nex$")

# remove ":" from names
names(binmat_nexus) <- str_replace_all(names(binmat_nexus), ":", ".")

# get as AAMultipleAlignment
binmat_biostrings <- 
  purrr::map(binmat_nexus, str_c, collapse = "") %>% 
  unlist() %>% 
  AAMultipleAlignment(.) 

# save in Phylip format which can be read by RAxML
write.phylip(x = binmat_biostrings, filepath = file.path(outpath, str_c(binmat_name, ".phy")))

