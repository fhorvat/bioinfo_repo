### INFO: 
### DATE: Fri Oct 18 12:52:33 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Analysis/piRNA_clusters")

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

# piRNA clusters found by PILFER
pirna_list <- list.files(inpath, pattern = ".*collapsed.clusters.txt", full.names = T)

######################################################## READ DATA
# read piRNA clusters table
pirna_tb <- purrr::map(pirna_list, function(path){
  
  # read table, get coordinates, add sample ID
  readr::read_delim(file = path, delim = "\t", col_names = c("coordinates", "score")) %>% 
    tidyr::separate(col = coordinates, into = c("seqnames", "coordinates"), sep = ":") %>% 
    tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-") %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.24to31nt\\.collapsed\\.clusters\\.txt"), 
                  genotype = str_extract(sample_id, "Mov10l_KO|Mov10l_WT"))
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
# get ranges
pirna_gr <- 
  pirna_tb %>% 
  GRanges(.) %>% 
  reduce(., ignore.strand = T) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(gene_id = str_c("pilfer_cl", 1:nrow(.)), 
                strand = "*") %>% 
  dplyr::select(GeneID = gene_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
  readr::write_delim(., path = file.path(outpath, "expression", "piRNA_clusters.testis_Mov10l.PILFER.saf"), delim = "\t")

  
 