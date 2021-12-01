#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: count reads over ENSEMBL
### DATE: Wed Sep 26 17:40:45 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/ensembl_counts/%EXPERIMENT")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genomes/%EXPERIMENT"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.93*.reducedExons.RDS$", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/%EXPERIMENT"

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# get list of bam files
bam_path <- list.files(path = mapped_path, pattern = "*.total.bam$", full.names = T)

######################################################## MAIN CODE
### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = FALSE,
                                           ignore.strand = TRUE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, str_c("ensembl.93.%EXPERIMENT.se.RDS")))
