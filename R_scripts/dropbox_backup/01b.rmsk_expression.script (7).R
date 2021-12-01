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

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# .bed path
rmsk_bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/whole_rmsk/expression"
rmsk_bed_path <- file.path(rmsk_bed_path, "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.bed")

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
bam_path <- args$bam_path
bam_name <- args$bam_name

# don't overwrite your .bam
if(bam_name == bam_path){
  stop("Something's wrong with your bam name, you don't want to overwrite your original bam file")
}

######################################################## READ DATA
# read .bed
rmsk_bed <- rtracklayer::import.bed(rmsk_bed_path)

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
         read_strand = as.character(strand(bam_hits)),
         read_width = width(bam_hits),
         rmsk_name = mcols(rmsk_hits)$name) %>%
  dplyr::group_by(read_name) %>%
  dplyr::sample_n(1) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(rmsk_name = str_remove(rmsk_name, ".*\\.")) %>%  
  dplyr::group_by(rmsk_name, read_width) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup(.)

# save as .RDS
saveRDS(overlaps_tb, file.path(outpath, "expression_sum.RDS_files", str_c(bam_name, ".expression.RDS")))
