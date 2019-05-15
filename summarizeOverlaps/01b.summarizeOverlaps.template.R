#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get expression of developmental profile
### DATE: /common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression/summarizedExperiments
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# set experiment name
experiment <- "%EXPERIMENT"

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# mapped path
mapped_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/", experiment, "/Data/Mapped/STAR_mm10")

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# get list of bam files
bam_path <- list.files(path = mapped_path, pattern = "*.total.bam$|*.genome.Aligned.sortedByCoord.out.bam$", full.names = T)

######################################################## MAIN CODE
### get pairing design of an experiment
# get bam names
bam_names <- 
  basename(bam_path) %>% 
  str_remove(., "\\.total\\.bam|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$")

# get pairing
pairing <- 
  str_extract(bam_names, "SE$|PE$") %>% 
  unique(.)

# sanity check
if(length(pairing) > 1 | !(pairing %in% c("SE", "PE"))) stop("Something's wrong with pairing of the data")

# check if single end
isSingleEnd <- ifelse(pairing == "SE", T, F)


### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = isSingleEnd,
                                           ignore.strand = TRUE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, basename(exons_path) %>% str_replace(., "reducedExons.RDS", str_c(experiment, ".se.RDS"))))
