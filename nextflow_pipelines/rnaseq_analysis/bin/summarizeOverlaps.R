#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get gene counts
### DATE: Mon Jun 17 15:26:08 2019
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats and mapped logs
bam_path <- args$mapped_bams_list
gtf_path <- args$gtf_path

######################################################## READ DATA
# read .gtf
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
### get pairing design of an experiment
# get bam names
bam_names <- 
  basename(bam_path) %>% 
  str_remove(., "\\.bam$")

# get pairing
pairing <- 
  str_extract(bam_names, "SE$|PE$") %>% 
  unique(.)

# sanity check
if(length(pairing) > 1 | !(pairing %in% c("SE", "PE"))) stop("Something's wrong with pairing of the data")

# check if single end
isSingleEnd <- ifelse(pairing == "SE", T, F)


### get reduced exons coordinates 
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(gtf_tb, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F)


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
saveRDS(se, file = file.path(outpath, basename(gtf_path) %>% str_replace(., "\\.gtf\\.gz$", ".se.RDS")))
