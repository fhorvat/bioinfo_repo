### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/rmsk")

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

library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LINE1s rmsk path
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.annotated_exons.csv")

# L1Base path
l1base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1Base/LINE1_annotation.GRanges.L1Base.RDS"

######################################################## READ DATA
# read LINE1 table
line1_tb <- readr::read_csv(line1_path)

# read L1Base
l1base_gr <- readRDS(l1base_path)

######################################################## MAIN CODE
### add info about L1Base annotation overlap
# to GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.) 


# find overlaps and percentage of overlap of each hit
line1_l1base_foverlaps <- findOverlaps(line1_gr, l1base_gr, ignore.strand = T)
# overlaps <- pintersect(line1_gr[queryHits(line1_l1base_foverlaps)], l1base_gr[subjectHits(line1_l1base_foverlaps)], ignore.strand = T)
# percentOverlap <- width(overlaps) / width(line1_gr[queryHits(line1_l1base_foverlaps)])

# extract ID's and genes
line1_l1base_tb <- 
  tibble(rmsk_id = line1_gr[queryHits(line1_l1base_foverlaps)]$rmsk_id,
         L1Base_id_overlap = l1base_gr[subjectHits(line1_l1base_foverlaps)]$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>%
  dplyr::summarise(L1Base_id_overlap = str_c(L1Base_id_overlap, collapse = ", "))

# add info about overlaping exons to table
line1_tb %<>% dplyr::left_join(., line1_l1base_tb, by = "rmsk_id")

# save table with info about exon overlap
line1_tb %T>%
  readr::write_csv(., file.path(outpath, "LINE1.4000nt_plus.ORFs.annotated_exons_and_L1Base.csv"))

