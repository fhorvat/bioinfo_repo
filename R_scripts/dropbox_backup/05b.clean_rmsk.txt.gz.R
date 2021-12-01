### INFO: reads rmsk.txt.gz, outputs tab-delim table with coordinates and repeat name/class/family
### DATE: 28 08. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/vfranke")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)

######################################################## SOURCE

######################################################## PATH VARIABLES
out_table <- list.files(path = getwd(), pattern = "rmsk\\..*\\.raw\\.fa\\.out\\.gz", full.names = T)

######################################################## MAIN CODE
sapply(out_table, function(x){
  
  ### clean and save repeatMasker
  readr::read_table2(file = x, skip = 3, col_names = F) %>% 
    dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_ID = X15) %>%
    tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
    dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
    readr::write_delim(str_replace(x, "raw", "clean"), delim = "\t")
  
  
})
