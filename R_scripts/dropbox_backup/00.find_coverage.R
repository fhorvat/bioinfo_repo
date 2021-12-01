### INFO: 
### DATE: Sat Oct 05 13:47:57 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression_coverage")

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

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Mapped/bbmap_mesAur1/4_merged_replicates"

# bigWig path
coverage_path <- 
  list.files(mapped_path, pattern = ".*\\.bw", full.names = T) %>% 
  .[!str_detect(., "scaled")]

######################################################## READ DATA
# read coverage from bigWig
coverage_list <- purrr::map(coverage_path, function(path){
  
  # read bigWig file
  coverage <- 
    rtracklayer::import(path) %>% 
    .[mcols(.)$score > 0] %>% 
    reduce(., ignore.strand = T) %>% 
    .[width(.) >= 100]
  
}) %>% 
  set_names(., coverage_path %>% basename(.) %>% str_remove_all(., "^s_GV_|\\.bw$"))

######################################################## MAIN CODE
# join, reduce
coverage_joined <- 
  purrr::reduce(coverage_list, c) %>% 
  GenomicRanges::reduce(., ignore.strand = T)

# filter table and save
coverage_tb <- 
  coverage_joined %>%
  as.data.table(.) %>% 
  .[, gene_id := str_c(seqnames, ":", start, "-", end)] %>%
  .[]

# save
coverage_tb %>%
  as_tibble(.) %>% 
  dplyr::select(gene_id, seqnames, start, end) %T>%
  readr::write_csv(., path = file.path(outpath, "joined_coverage.perfect_multimappers.coordinates.csv"))

# save as SAF
coverage_tb %>% 
  as_tibble(.) %>% 
  dplyr::select(GeneID = gene_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
  readr::write_delim(., path = file.path(outpath, "joined_coverage.perfect_multimappers.coordinates.saf"), delim = "\t")
