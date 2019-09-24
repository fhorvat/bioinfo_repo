#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get expression of developmental profile
### DATE: /common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression/summarizedExperiments
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
ensembl_version <- args$ensembl_version
genome_path <- args$genome_path
mapped_path <- args$mapped_path
threads <- args$threads


### other experiment paths
# set base experiment path
base_path <- file.path(mapped_path, "../..")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression")


### bam files
# get list of bam files
bam_path <- list.files(path = mapped_path, pattern = "*.bam$", full.names = T)


### genome
# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
### get pairing design of an experiment
# get bam names
bam_names <-
  basename(bam_path) %>%
  str_remove(., "\\.total\\.bam|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$|\\.bam$")

# get pairing
pairing <-
  str_extract(bam_names, "SE$|PE$") %>%
  unique(.)

# sanity check
if(length(pairing) > 1 | !(pairing %in% c("SE", "PE"))) stop("Something's wrong with pairing of the data")

# check if single end
isSingleEnd <- ifelse(pairing == "SE", T, F)


### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = threads))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = isSingleEnd,
                                           ignore.strand = TRUE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, basename(exons_path) %>% str_replace(., "reducedExons.RDS", str_c(experiment, ".se.RDS"))))

