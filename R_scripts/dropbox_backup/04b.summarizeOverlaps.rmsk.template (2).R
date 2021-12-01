#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarizes counts in developmental profile over LTRs and LINE1s from repeatMasker
### DATE: Mon Feb 25 14:20:36 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/summarizedExperiments.all_LINEs_LTRs")

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
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeatMasker VIZ path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

# accessory datasets path
accessory_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq"

# experiment path
experiment_path <- file.path(accessory_path, "%EXPERIMENT", "Data/Mapped/STAR_mm10")

# bam paths
bam_path <- list.files(experiment_path, pattern = ".*\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$", full.names = T)

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
### prepare data
# form GRanges
rmsk_gr <- 
  rmsk_tb %>% 
  dplyr::filter(repClass %in% c("LINE", "LTR")) %>% 
  GRanges(.)
names(rmsk_gr) <- rmsk_gr$rmsk_id

### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps - single end
se <-
  GenomicAlignments::summarizeOverlaps(features = rmsk_gr,
                                       reads = bamfiles,
                                       param = ScanBamParam(tagFilter = list("nM" = 0)),
                                       mode = "Union",
                                       singleEnd = TRUE,
                                       ignore.strand = TRUE) %>%
  saveRDS(., file = file.path(outpath, str_c("rmsk.L1_and_LTRs.all.", "%EXPERIMENT", ".perfect.20190225.se.RDS")))
