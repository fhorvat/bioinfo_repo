### INFO: 
### DATE: Thu Oct 17 14:13:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
bed_path <- args$bed_path

# bed name
bed_name <- bed_path %>% basename(.) %>% str_remove(., "\\.bed$")

######################################################## READ DATA

######################################################## MAIN CODE
# read, collapse reads, save as bed with count of alignments as score
data.table::fread(bed_path, col.names = c("seqnames", "start", "end", "qname", "score", "strand")) %>% 
  .[, .(score = .N, qname = qname[1]), by = c("seqnames", "start", "end", "strand")] %>% 
  setcolorder(., c("seqnames", "start", "end", "qname", "score", "strand")) %T>% 
  readr::write_delim(., file.path(outpath, str_c(bed_name, ".collapsed.bed")), delim = "\t", col_names = FALSE)
