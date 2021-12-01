### INFO: 
### DATE: Tue Sep 24 20:46:38 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1")

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
inpath <- getwd()

# set outpath
outpath <- getwd()

# repeatMasker path
rmsk_path <- list.files(path = inpath, pattern = "rmsk\\..*\\.raw\\.fa\\.out\\.gz", full.names = T)

# scaffold names path
scaffolds_path <- "/common/DB/genome_reference/golden_hamster/vfranke/scaffold_names.txt"

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

# read scaffold names
scaffolds_tb <- readr::read_delim(scaffolds_path, delim = "\t")

######################################################## MAIN CODE
# clean rmsk, join with other scaffolds names, save clean
rmsk_tidy <- 
  rmsk_tb %>% 
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_ID = X15) %>%
  dplyr::left_join(., scaffolds_tb %>% dplyr::select(`RefSeq Accession.version`, `GenBank Accession.version`), by = c("seqnames" = "RefSeq Accession.version")) %>% 
  dplyr::select(-seqnames) %>% 
  dplyr::select(seqnames = `GenBank Accession.version`, everything()) %>% 
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
  readr::write_delim(., str_replace(rmsk_path, "raw", "clean"), delim = "\t")

# join raw with other scaffold names, save
rmsk_raw <- 
  rmsk_tb %>% 
  dplyr::left_join(., scaffolds_tb %>% dplyr::select(`RefSeq Accession.version`, `GenBank Accession.version`), by = c("X5" = "RefSeq Accession.version")) %>% 
  dplyr::select(-X5) %>% 
  dplyr::select(X1:X4, X5 = `GenBank Accession.version`, everything()) %T>% 
  readr::write_delim(., str_replace(rmsk_path, "raw", "RefSeq.raw"), delim = "\t", col_names = F)

  

