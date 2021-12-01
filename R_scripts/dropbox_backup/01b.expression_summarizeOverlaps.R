#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
# set experiment
experiment <- "hamster_oocyte_Mov10l.RNAseq"
single_end <- TRUE

# set working dir
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/FPKM_tables/L1_joined_rmskID", experiment))

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

# documentation path 
documentation_path <- file.path(inpath, "../../../Documentation/L1_joined_rmskID/mesAur1.L1_joined_rmskID_20190929.FH")

# list of LINE1 full length elements path
line1_coords_path <- file.path(documentation_path, "LINE1_whole.joined_rmsk_ID.csv")

# datasets path
datasets_path <- file.path(inpath, "../../../datasets", experiment, "Mapped/mesAur1.L1_joined_rmskID_20190929.FH")

# get list of bam files
bam_path <- list.files(path = datasets_path, pattern = "*.bam$", full.names = T)

######################################################## READ DATA
# read rmsk list of LINE1 full length elements
line1_tb <- readr::read_csv(line1_coords_path)

######################################################## MAIN CODE
### prepare data
# form GRanges
line1_gr <- GRanges(line1_tb)
names(line1_gr) <- line1_gr$rmsk_id

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
  saveRDS(., file = file.path(outpath, str_c("LINE1_whole.joined_rmsk_ID.", experiment, ".se.RDS")))

