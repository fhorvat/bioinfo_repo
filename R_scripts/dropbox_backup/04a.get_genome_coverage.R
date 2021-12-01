### INFO: 
### DATE: Sat Oct 05 13:47:57 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/miR-205_pig/genomic_targets")

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

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/pig.susScr11/4_merged_replicates"

# bigWig path
coverage_path <- file.path(mapped_path, "s_Pig_all.PE.bw")

######################################################## READ DATA
# read coverage from bigWig
coverage <- rtracklayer::import(coverage_path)

######################################################## MAIN CODE
# filter table and save
coverage_tb <- 
  coverage %>%
  .[mcols(.)$score > 0] %>% 
  reduce(., ignore.strand = T) %>% 
  as.data.table(.) %>% 
  .[, gene_id := str_c(seqnames, ":", start, "-", end)] %>%
  .[]

# save
coverage_tb %>%
  as_tibble(.) %>% 
  dplyr::select(gene_id, seqnames, start, end) %T>%
  readr::write_csv(., path = file.path(outpath, "s_Pig_all.PE.joined_coverage.coordinates.csv"))

# # save as SAF
# coverage_tb %>% 
#   as_tibble(.) %>% 
#   dplyr::select(GeneID = gene_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
#   readr::write_delim(., path = file.path(outpath, "s_Pig_all.PE.joined_coverage.coordinates.saf"), delim = "\t")

# export as .gtf
coverage_gr <- GRanges(coverage_tb)
names(coverage_gr) <- mcols(coverage_gr)$gene_id
export.gff3(object = coverage_gr, con = file.path(outpath, "s_Pig_all.PE.joined_coverage.coordinates.gff"))


