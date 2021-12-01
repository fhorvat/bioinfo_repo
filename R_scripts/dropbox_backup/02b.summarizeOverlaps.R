#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads over features (LINE1 elements)
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working dir
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
library(openxlsx)
library(purrr)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### IN AND OUT
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### COMMAND LINE ARGUMENTS
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)
dataset_path <- args$dataset_path
feature_coordinates <- args$feature_coordinates


### DATASET
# bam path
bam_path <- list.files(path = dataset_path, pattern = "*.bam$", full.names = T)


######################################################## READ DATA
# read LINE1 coordinates
line1_gr <- readRDS(feature_coordinates)

######################################################## MAIN CODE
### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = threads))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = line1_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = single_end,
                                           ignore.strand = TRUE,
                                           fragments = !single_end) %>%
  saveRDS(., file.path(outpath, str_c((feature_coordinates %>% basename(.) %>% str_remove(., "\\.GRanges\\.RDS$")), 
                                      experiment, 
                                      "se.RDS", 
                                      sep = ".")))
