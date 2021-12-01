#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# experiment
experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation"

# LINE1s coordinates
line1_path <- file.path(documentation_path, "LINE1.4000nt_plus.ORFs.annotated_exons.csv")

# basepath
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression"

# datasets path
datasets_path <- file.path(base_path, "datasets", experiment, "Mapped/perfect_alignments.all_multimappers")

# get list of bam files
bam_path <- list.files(path = datasets_path, pattern = "*.bam$", full.names = T)

######################################################## READ DATA
# read LINE1 table
line1_tb <- 
  readr::read_csv(line1_path) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

######################################################## MAIN CODE
### prepare data
# create GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.)
names(line1_gr) <- mcols(line1_gr)$rmsk_id

### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = threads))

# summarize overlaps - single end
se <- GenomicAlignments::summarizeOverlaps(features = line1_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = single_end,
                                           ignore.strand = TRUE,
                                           fragments = !single_end) %>%
  saveRDS(., file = file.path(outpath, str_c(basename(line1_path) %>% str_remove(., "\\.csv$"), experiment, "se.RDS", sep = ".")))
