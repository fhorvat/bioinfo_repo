### INFO:
### DATE: Thu Jul 30 16:34:41 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

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
library(purrr)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

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

# arguments from command line
bam_path <- args$bam_path
bam_name <- args$bam_name
bed_path <- args$bed_path

# don't overwrite your .bam
if(bam_name == bam_path){
  stop("Something's wrong with your bam name, you don't want to overwrite your original bam file")
}

######################################################## READ DATA
# read .bed
rmsk_bed <- rtracklayer::import.bed(bed_path)

# read alignments from .bam
bam_alignments <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T)

######################################################## MAIN CODE
### for each alignment find one random sense and antisense hit
# find overlaps with repeatMasker for each alignment
bam_gr <-
  bam_alignments %>%
  unlist(.)

# find all overlaps
overlaps <- findOverlaps(bam_gr, rmsk_bed, ignore.strand = T)

# get hits
bam_hits <- bam_gr[queryHits(overlaps)]
rmsk_hits <- rmsk_bed[subjectHits(overlaps)]

# summarize expression per repName, read width and sense
overlaps_tb <-
  tibble(read_name = names(bam_hits),
         rmsk_name = mcols(rmsk_hits)$name) %>%
  dplyr::group_by(read_name) %>%
  dplyr::sample_n(1) %>%
  dplyr::ungroup(.) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup(.)

# save as .RDS
saveRDS(overlaps_tb, file.path(outpath, "expression_sum.RDS_files", str_c(bam_name, ".expression.RDS")))

