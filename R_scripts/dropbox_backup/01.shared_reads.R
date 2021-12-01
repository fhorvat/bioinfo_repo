### INFO: 
### DATE: Thu Oct 10 04:42:16 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/family_sequences")

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

# LINE1 coordinates path
bed_path <- file.path(inpath, "L1Md_all.bed")

# filtered bam file path
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Dicer_Mili_KO/Mapped/perfect_alignments.all_multimappers/4_merged_replicates/L1_families_subset/s_GV_DBL.L1Md.bam"

######################################################## READ DATA
# read bam
bam <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T, param = ScanBamParam(what = c("qname")))
  
######################################################## MAIN CODE
