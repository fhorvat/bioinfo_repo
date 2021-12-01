### INFO: count small RNA reads mapped to antisense features in repeatMasker 
### DATE: Tue Jul 31 23:18:59 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/repeat_expression.20180730/antisense_counts_rmsk")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# wideScreen()

######################################################## FUNCTIONS
# class alignments hierarchically 
countAntisenseReads <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA){
  
  # initialize count vector
  count_vector <- rep(0, length(subject_ranges))

  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <- GenomicRanges::grglist(chunk) 
    
    # count overlaps between two GRanges - chunk of alignments and annotation
    hits_all <- countOverlaps(subject_ranges, chunk, ignore.strand = T)
    
    # count sense overlaps
    hits_sense <- countOverlaps(subject_ranges, chunk, ignore.strand = F)
    
    # calculate antisense overlaps
    hits_antisense <- hits_all - hits_sense
    
    # sum chunk of antisense overlaps with count vector
    count_vector <- count_vector + hits_antisense
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return vector with counts
  return(count_vector)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get bam file path, name, experiment
# bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Tam_2008_Nature_GSE10364/Data/Mapped/STAR_mm10_new/s_oocyte_19to30.SE.genome.Aligned.sortedByCoord.out.bam"
bam_path <- "%BAM_PATH"
bam_name <- basename(bam_path) %>% str_remove_all(., ".SE.genome.Aligned.sortedByCoord.out.bam|.SE|.bam")
experiment_name <- dirname(bam_path) %>% str_remove_all("/common.*small_RNAseq/|/Data.*$")

######################################################## READ DATA
# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
# class antisense reads over repeatMasker
antisense_count <- countAntisenseReads(bam_path = bam_path, 
                                       subject_ranges = rmsk_gr, 
                                       yield = 1000000, 
                                       isFirstInPair = NA)

# save as RDS
saveRDS(object = antisense_count, file = file.path(outpath, str_c("antisense_count.rmsk.", experiment_name, bam_name, "RDS", sep = ".")))
