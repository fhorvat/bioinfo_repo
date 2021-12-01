#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/FPKM_tables/Dicer_Mili_KO")

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

# list of LINE1 full length elements path
line1_coords_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

### 
# set experiment
experiment <- "Dicer_Mili_KO"
single_end <- FALSE

# datasets path
datasets_path <- file.path(inpath, "../../datasets", experiment, "Mapped/mm10_masked")

# get list of bam files
bam_path <-
  list.files(path = datasets_path, pattern = "*.bam$", full.names = T) %>%
  .[str_detect(., "mm10_masked")] %>%
  .[!str_detect(., "merged|almost_perfect")]

######################################################## READ DATA
# read Zoe's list of LINE1 full length elements
line1_coords <- read_csv(line1_coords_path)

######################################################## MAIN CODE
### prepare data
# form GRanges
line1_gr <- GRanges(line1_coords)
names(line1_gr) <- line1_gr$id

### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))

# summarize overlaps - single end
se <- GenomicAlignments::summarizeOverlaps(features = line1_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = single_end,
                                           ignore.strand = TRUE) %>%
  saveRDS(., file = file.path(outpath, str_c("L1s_nested_ours_20180516.ZJM.tidy.", experiment, ".se.RDS")))
