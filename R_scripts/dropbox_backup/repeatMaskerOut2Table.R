#!/home/students/fhorvat/R/bin/Rscript
### qsub -V -q MASTER -l select=1:mem=5gb -N rmskOutToTable -j oe -o /common/WORK/fhorvat/reference/cow/bosTau8 /common/WORK/fhorvat/reference/cow/bosTau8/repeatMaskerOut2Table.R
### INFO: reads .out.gz from repeatMasker, outputs tab-delim table with coordinates and repeat name/class/family
### DATE: 28 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/reference/cow/bosTau8")

######################################################## PATH VARIABLES
out_table <- list.files(path = getwd(), pattern = "*.fa.out.gz", full.names = T)
table_name <- 
  str_replace_all(out_table, "\\/.*\\/|.fa.out.gz", "") %>% 
  str_c(., ".rmsk.txt.gz")
  
######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)

######################################################## MAIN CODE
readr::read_table(file = out_table, skip = 3, col_names = F) %>% 
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11) %>% 
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>% 
  dplyr::mutate(strand = replace(strand, strand == "C", "*")) %T>%
  readr::write_delim(file.path(getwd(), table_name))
