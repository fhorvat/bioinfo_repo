### INFO: 
### DATE: Wed Jul 14 16:00:28 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/Mollusca/Arion_vulgaris.Schrodl")

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

# rmsk gff path
rmsk_gff_path <- file.path(inpath, "Arion.TEs.gff")

######################################################## READ DATA
# read repeatMasker gff
rmsk_gff <- rtracklayer::import(rmsk_gff_path)

######################################################## MAIN CODE
# seqnames, start, end, strand, repName, repClass, repFamily, rmsk_id
rmsk_tb <- 
  as_tibble(rmsk_gff) %>% 
  dplyr::mutate(repName = str_remove_all(Target, "lcl\\|| .*")) %>% 
  tidyr::separate(col = Class, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(rmsk_id = 1:n()) %>% 
  dplyr::select(seqnames, start, end, strand, repName, repClass, repFamily, rmsk_id)

# save
readr::write_delim(rmsk_tb, file.path(outpath, "rmsk.AriVul.20210714.clean.fa.out.gz"), delim = "\t")

  
