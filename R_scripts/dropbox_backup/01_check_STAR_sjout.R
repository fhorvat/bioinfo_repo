#!/home/students/fhorvat/R/bin/Rscript
### INFO: checks genomic position for overlap with annotated and non-canonical splice junctions from STAR output
### DATE: 14. 09. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/mouse/Data/Mapped/STAR_mm10/SINE_region")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## PATH VARIABLES
outpath <- getwd()
sj_list <- list.files(path = outpath, pattern = "*SJ.out.tab$", full.names = T)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "parseSJOut.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# read SJ.out table from STAR
sjout <-
  parseSJOut(sj_list[3]) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

######################################################## MAIN CODE
# findOverlaps with feature of interest chr19:5801592-5802869
SINE_gr <- GenomicRanges::GRanges(seqnames = "chr19", ranges = IRanges::IRanges(start = 5801592, end = 5802869), strand = "*")
SINE_overlaps <- findOverlaps(SINE_gr, sjout, ignore.strand = T)
sjout[subjectHits(SINE_overlaps)]


