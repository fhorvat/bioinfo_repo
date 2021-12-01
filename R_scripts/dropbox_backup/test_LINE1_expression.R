#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 04. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# counts reads over ranges
countReads <- function(granges_bam, pair_end = pair_end_seq, query_range){
  
  # overlap reads
  overlap_reads <- IRanges::subsetByOverlaps(x = granges_bam, ranges = query_range, type = "within", ignore.strand = T)
  
  # add column to sum table
  reads_sum <- tibble(n_read = length(unique(names(overlap_reads))), 
                      n_alignment = length(overlap_reads))
  
  return(reads_sum)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get path of bam file
bam_path <- "%BAM_PATH"

# get name of the samples
bam_name <-
  basename(bam_original_path) %>%
  stringr::str_replace_all(., ".genome.Aligned.sortedByCoord.out.bam", "")

# determine whether .bam comes from pair-end sequencing
pair_end_seq <- str_detect(string = bam_name, pattern = "PE")

######################################################## READ DATA
# coordinates
coordinates_gr <- readRDS(file = coordinates_path)

######################################################## MAIN CODE
### count reads over features hierarchically - complete match: rDNA -> repeat -> exon -> other
# read pair-end or single end bam file, if pair-end set unique names for each mate in fragment 
if(pair_end_seq){
  
  # read alignment, set unique names for read 1 and 2 in pair, unlist to GRanges
  gbam_original <- bamToGRangesList(bam_original_path) 
  names(gbam_original) <- str_c(names(gbam_original), 1:2)
  gbam_original <- unlist(gbam_original)
  
  # get counts of reads and alignments mapped to coordinate
  gbam_counts <- classReads(granges_bam = .)
  
}else{
  
  # read alignment, unlist to GRanges
  gbam_original <-
    bamToGRangesList(bam_original_path) %>%
    unlist(.)
  
  # get counts of reads and alignments mapped to coordinate
  gbam_counts <- classReads(granges_bam = .)
  
}


