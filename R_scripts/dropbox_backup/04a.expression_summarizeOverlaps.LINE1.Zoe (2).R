### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression")

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
line1_coords_path <- file.path(inpath, "Documentation", "LINE1_full_length.20180517.ZJM.tidy.csv")

# datasets path
datasets_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/datasets" 

# get list of bam files
bam_path <- 
  list.files(path = datasets_path, pattern = "*.bam$", full.names = T, recursive = T) %>% 
  .[str_detect(., "mm10_masked")] %>% 
  .[!str_detect(., "PWD_CAST")]

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
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps - single end
se_SE <- 
  GenomicAlignments::summarizeOverlaps(features = line1_gr,
                                       reads = bamfiles[str_detect(names(bamfiles), "SE")],
                                       mode = "Union",
                                       singleEnd = TRUE,
                                       ignore.strand = TRUE) %>% 
  saveRDS(., file = file.path(outpath, "results", "LINE1_full_length.20180517.ZJM.tidy.SE.se.RDS"))

# summarize overlaps - paired end
se_PE <- 
  GenomicAlignments::summarizeOverlaps(features = line1_gr,
                                       reads = bamfiles[str_detect(names(bamfiles), "PE")],
                                       mode = "Union",
                                       singleEnd = FALSE,
                                       ignore.strand = TRUE) %>% 
  saveRDS(., file = file.path(outpath, "results", "LINE1_full_length.20180517.ZJM.tidy.PE.se.RDS"))



