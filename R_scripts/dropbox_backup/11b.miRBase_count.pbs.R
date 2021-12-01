### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Wed Nov 28 15:03:32 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRbase path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# get bam file path, name, experiment
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR.new"

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*bam$", full.names = T)

######################################################## READ DATA
# read miRbase gtf
mirbase_gr <- rtracklayer::import.gff(con = mirbase_path) 

######################################################## MAIN CODE
# get ranges of mature miRNA
mirna_gr <- mirbase_gr[mcols(mirbase_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name", "Derives_from")]

### get counts of reads over miRNA
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_paths, yieldSize = 3000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = mirna_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = FALSE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, "miRbase.Eliska_mESC_MosIR.se.RDS"))


