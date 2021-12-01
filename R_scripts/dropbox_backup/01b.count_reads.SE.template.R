### INFO: Count reads and get FPKM of genes in mESC and oocytes sequenced in February 2018 and in June 2018
### DATE: Fri Nov 30 14:06:13 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Analysis/expression")

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
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation/mESC_oocytes_2018.sample_table.csv"

# reduced exons path
exons_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.reducedExons.RDS"

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# filter sample table
sample_table %<>%
  dplyr::filter(str_detect(sample_id, "^s_ESC_DX_.*"))

# get count of reads, save summarizedExperiment as RDS
bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = TRUE)

saveRDS(se, file = file.path(outpath, "mESC_DX.GRCm38.91.reducedExons.summarizedOverlaps.RDS"))
