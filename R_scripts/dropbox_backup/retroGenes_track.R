#!/home/students/fhorvat/R/bin/Rscript
### INFO: read GFF file, for each gene take transcipt with most exons
### DATE: 22. 5. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(readr)
library(magrittr)
library(dplyr)
library(stringr)
library(tibble)
library(GenomicRanges)

######################################################## PATH VARIABLES
table_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_RetroGenesInfoV6_20170602.txt.gz"

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
# read .gtf
retroGenes <- 
  read_delim(file = table_path, delim = "\t", col_names = T, col_types = cols(.default = "c")) %>% 
  dplyr::select(seqnames = 1, start = chromStart, end = chromEnd, strand, name, type)
  
######################################################## MAIN CODE


