#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: finds coverage of repetitive elements from list and outputs total coverage ratio and longest coverage stretch in each element
### DATE: Sat Oct 05 13:47:57 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working dir
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

bw_path='../bw_subset/s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE.bw'
bw_name='s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE '
bed_path='../../../IAP.FLI_elements.bed'

######################################################## READ DATA
# read .bed
rmsk_bed <- rtracklayer::import.bed(bed_path)

# read bigWig file
element_coverage_full <- rtracklayer::import(bw_path)

######################################################## MAIN CODE
### prepare files
# create element table
rmsk_tb <- 
  rmsk_bed %>% 
  as_tibble(.)


### find coverage for each element
# reduce coverage
element_coverage_full <- reduce(element_coverage_full, ignore.strand = T)

# find overlaps between two GRanges
hits <- findOverlaps(rmsk_bed, element_coverage_full, ignore.strand = T)

# extract all overalaping features from subject as list
element_coverage <- extractList(element_coverage_full, as(hits, "List"))

# intersect with LINE1 coordinates
element_coverage <- pintersect(element_coverage, rmsk_bed, ignore.strand = T)

# set names
names(element_coverage) <- mcols(rmsk_bed)$name 


### find longest stretch and percentage of coverage for each element
# unlist
element_coverage_tb <- 
  element_coverage %>% 
  unlist(.)
mcols(element_coverage_tb)$name <- names(element_coverage_tb)

# convert to table, summarize 
element_coverage_sum <- 
  element_coverage_tb %>% 
  as_tibble(.) %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(coverage_total = sum(width), 
                   coverage_longest = max(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., rmsk_tb %>% dplyr::select(name, width), by = "name") %>% 
  dplyr::mutate(total_coverage_ratio = round((coverage_total / width), 3), 
                longest_coverage_ratio = round((coverage_longest / width), 3)) %>% 
  dplyr::select(name, total_coverage_ratio, longest_coverage_ratio) %>% 

# join to one table
line1_coverage_tb <- 
  purrr::reduce(c(list(line1_tb %>% dplyr::select(rmsk_id)), line1_coverage_list), left_join, by = "rmsk_id") %>% 
  dplyr::mutate_at(vars(contains("coverage_ratio")), list(~ifelse(is.na(.), 0, .))) %>% 
  dplyr::mutate_at(vars(contains("longest_")), list(~ifelse(is.na(.), 0, .))) %T>%
  readr::write_csv(., path = file.path(outpath, str_c((feature_coordinates %>% basename(.) %>% str_remove(., "\\.GRanges\\.RDS$")), 
                                                      experiment, 
                                                      "coverage.csv", 
                                                      sep = ".")))
