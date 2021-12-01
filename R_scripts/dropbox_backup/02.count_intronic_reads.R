### INFO: 
### DATE: Fri Jun 28 16:41:11 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Fugaku_intron_CpG/count_intronic_reads")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "../intronic_reads")

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# intronic .bed path
introns_bed_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.bed")

######################################################## READ DATA
# read intron coordinates from .bed
introns_gr <- rtracklayer::import.bed(introns_bed_path)

# get .bam path
bam_path <- list.files(inpath, pattern = ".intronic.bam", full.names = T)

######################################################## MAIN CODE
# get unique names of introns
mcols(introns_gr)$score <- NULL
mcols(introns_gr)$name <- make.unique(mcols(introns_gr)$name)
names(introns_gr) <- mcols(introns_gr)$name
mcols(introns_gr)$name <- NULL

### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = introns_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = F,
                                           ignore.strand = TRUE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.se.RDS"))



