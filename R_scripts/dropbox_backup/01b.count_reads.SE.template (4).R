### INFO: get expression of all genes in lncKO data
### DATE: Wed May 23 19:20:31 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/diffExp")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(tibble)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/diffExp/lncKO.all.sample_table.csv"

######################################################## READ DATA
# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# get count of reads, save summarizedExperiment as RDS
bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
names(bamfiles) <- sample_table$sample_id
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = TRUE)

# save RDS
saveRDS(se, file = file.path(outpath, "ensembl.91.GRCm38.p5.20180512.lncKO.summarizedOveralaps.RDS"))
