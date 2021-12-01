### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/05_overlap_with_5pUTR_LINE1s")

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
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# minimap results path
bam_path <- file.path(inpath, "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.no_gaps.bam")

######################################################## READ DATA
# read mapping results
bam_gr <-  GenomicAlignments::readGAlignmentsList(file = bam_path, param = ScanBamParam(tag = "NM"), use.names = TRUE)

######################################################## MAIN CODE
# get ranges for each alignment
line1_tb <- 
  bam_gr %>% 
  unlist(.) 

# set read name as query 
mcols(line1_tb)$query_id <- names(line1_tb)
  
# get the alignment without mismatches for each hit
line1_tb %<>% 
  as_tibble(.) %>% 
  dplyr::filter(NM == 0) 

# save as table
line1_tb %>% 
  dplyr::select(seqnames, start, end, strand, query_id) %>% 
  readr::write_csv(., "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.coordinates.csv")

