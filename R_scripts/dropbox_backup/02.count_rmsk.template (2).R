#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts in developmental profile over LTRs and LINE1s from repeatMasker
### DATE: Thu Mar 07 21:58:57 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# bam paths
bam_path <- "%BAM_PATH"
bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")

# repeatMasker coordinates path
rmsk_coords_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/rmsk.L1_and_LTRs.filtered.20190306.removed_genes.GRanges.RDS"

######################################################## READ DATA
# read repeatMasker GRanges
rmsk_gr <- readRDS(rmsk_coords_path)

# read bam
bam_gr <- readGAlignmentsList(bam_path, param = ScanBamParam(what = "qname"))

######################################################## MAIN CODE
# split GAlignments list to aligments of same read
bam_gr_split <- 
  bam_gr %>% 
  unlist(.) %>% 
  split(., mcols(.)$qname)

# split GRanges to GRangesList
rmsk_gr_list <- split(rmsk_gr, rmsk_gr$class)

# count reads overlaping with repeatMasker = each read is counted only once for each of the classes in repeatMasker 
overlaps_table <- 
  findOverlaps(bam_gr_split, rmsk_gr_list, ignore.strand = T, type = "within") %>% 
  as.data.table(.) %>% 
  .[, list(count = .N), by = subjectHits] %>% 
  as.tibble(.) %>% 
  dplyr::mutate(rmsk_id = names(rmsk_gr_list)[subjectHits], 
                sample_id = bam_name) %>% 
  dplyr::select(rmsk_id, count, sample_id) %T>%
  saveRDS(., file = file.path(outpath, str_c("rmsk_counts.", bam_name, ".RDS")))
  