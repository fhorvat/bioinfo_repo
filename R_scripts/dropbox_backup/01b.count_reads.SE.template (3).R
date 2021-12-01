### INFO: Count reads and get FPKM of genes
### DATE: Tue Dec 11 22:20:10 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression/joined_with_2017")

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
# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.2017Sep_2018Dec.sampleTable.clean.csv"

# reduced exons path
exons_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.reducedExons.RDS"

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# get count of reads, save summarizedExperiment as RDS
bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam(workers = 13))
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = TRUE)

# save
saveRDS(se, file = file.path(outpath, "Lnc1_KO.2018Dec_2017Sep.GRCm38.91.reducedExons.summarizedOverlaps.RDS"))
