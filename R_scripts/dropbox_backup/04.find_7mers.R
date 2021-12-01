### INFO: 
### DATE: Mon Oct 29 16:48:51 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# all 3'UTR sequences 
aligned_3UTRs_path <- file.path(accessory_data_path, "aligned_UTR_sequences.mm.rn.bt.hs.RDS")

# 7mers path
seed_7mers_path <- file.path(accessory_data_path, "seeds.7mers_patterns.RDS")

# binned expression path
exp_bin_path <- file.path(inpath, "grid.Su.mas5.oocyte.50.bins.csv")

######################################################## READ DATA
# read aligned 3'UTRs sequences in mouse, rat, cow and human
aligned_3UTRs <- readRDS(aligned_3UTRs_path)

# read 7mers table
seed_7mers <- readRDS(seed_7mers_path)

# read binned expression
exp_bin <- readr::read_csv(exp_bin_path)

######################################################## MAIN CODE
# get mouse UTRs
mouse_UTRs <- 
  aligned_3UTRs[["Mus_musculus"]] %>% 
  DECIPHER::RemoveGaps(.) %>% 
  as.character(.) %>% 
  magrittr::set_names(., str_remove(names(.), "\\..*|\\|Mus_musculus")) %>% 
  .[names(.) %in% exp_bin$gene_id]
  
# join UTRs in the same relative expression bins
binned_UTRs <- purrr::map(1:50, function(bin){
  
  # genes in bin
  bin_genes <- 
    exp_bin %>% 
    filter(bin_relative == bin) %$%
    gene_id
  
  # join UTRs of binned genes
  bin_UTR <- 
    mouse_UTRs[names(mouse_UTRs) %in% bin_genes] %>% 
    str_c(., collapse = "N")
  
})

# get septamer patterns as a list
septamer_patterns <- list(`7mer-m8` = seed_7mers$seed_7mer_m8, 
                          `7mer-1a` = seed_7mers$seed_7mer_1a)

# count septamers in binned UTRs
binned_UTRs_m8_matches <- purrr::map(septamer_patterns$`7mer-m8`, function(septamer) stringr::str_count(string = binned_UTRs, pattern = septamer)) 

# count septamers in binned UTRs
binned_UTRs_m8_matches <- purrr::map(binned_UTRs, function(utr) stringr::str_count(string = utr, pattern = septamer_patterns$`7mer-m8`)) 


# count matches



