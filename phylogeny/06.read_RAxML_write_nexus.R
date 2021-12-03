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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list files
tree_list_path <- list.files(inpath, "RAxML_fastTreeSH_Support.Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences_PS_CDD..*.msa.fasta.shsupport")

######################################################## READ DATA
# read files
tree_list <- purrr::map(tree_list_path, function(path){
  
  # out name
  tree_name <- 
    path %>% 
    basename(.) %>% 
    str_replace(., "\\.shsupport$", ".tre")
  
  # read as tree and save as nexus
  ape::read.tree(file = path) %>% 
    ape::write.nexus(., file = file.path(outpath, tree_name))
  
})

######################################################## MAIN CODE

