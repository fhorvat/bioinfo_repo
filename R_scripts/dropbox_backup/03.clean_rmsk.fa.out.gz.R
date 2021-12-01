### INFO: reads .out.gz from repeatMasker, outputs tab-delim table with coordinates and repeat name/class/family
### DATE: 28 08. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(purrr)

######################################################## SOURCE

######################################################## PATH VARIABLES
out_table <- list.files(path = getwd(), pattern = "rmsk.*raw.fa.out.gz", full.names = T, recursive = T)

######################################################## MAIN CODE
purrr::map(out_table, function(path){
  
  cat("\n", path, "\n", "\n")
  
  readr::read_table2(file = path, skip = 3, col_names = F) %>%
    dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_ID = X15) %>%
    tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
    dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
    readr::write_delim(str_replace(path, "raw", "clean"), delim = "\t")
  
})
