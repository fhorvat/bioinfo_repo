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
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/oocyte_WT/CNOT6L/Mapped/mm10_masked.perfect/4_merged_tracks/s_GV_WT.PE.perfect.merged.bam"
bam_name <- basename(bam_path) %>% str_remove(., "\\.perfect\\.merged\\.bam")
unmasked_bam_path <- file.path(dirname(bam_path), "filtered_reads", str_c(bam_name, ".unmasked.filtered.bam"))

######################################################## READ DATA
# read bam
bam_masked <- readGAlignmentsList(bam_path, param = ScanBamParam(what = "qname", tag = "AS"))
bam_unmasked <- readGAlignmentsList(unmasked_bam_path, param = ScanBamParam(what = "qname", tag = "AS"))

######################################################## MAIN CODE
### get maximum alignment score for each read
# masked bam
bam_masked_tb <- 
  as.tibble(bam_masked) %>% 
  group_by(group) %>% 
  summarize(qname = unique(qname), 
            AS = sum(AS)) %>% 
  group_by(qname) %>% 
  summarize(AS_masked = max(AS)) %>% 
  ungroup(.)

# unmasked bam
bam_unmasked_tb <- 
  as.tibble(bam_unmasked) %>% 
  group_by(group) %>% 
  summarize(qname = unique(qname), 
            AS = sum(AS)) %>% 
  group_by(qname) %>% 
  summarize(AS_unmasked = max(AS)) %>% 
  ungroup(.)

# join together, get read names which have higher alignment score in unmasked mapping 
alignment_scores <- 
  left_join(bam_masked_tb, bam_unmasked_tb, by = "qname") %>% 
  
  

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
