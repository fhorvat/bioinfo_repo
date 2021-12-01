#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarizes counts in developmental profile over LTRs and LINE1s from repeatMasker
### DATE: Wed Mar 06 11:05:00 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/summarizedExperiments.all_LINEs_LTRs")

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

# experiment name 
experiment_name <- "%EXPERIMENT"

# accessory datasets path
accessory_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq"

# experiment path
experiment_path <- file.path(accessory_path, experiment_name, "Data/Mapped/STAR_mm10")

# bam paths
bam_path <- list.files(experiment_path, pattern = ".*\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$", full.names = T)

# repeatMasker coordinates path
rmsk_coords_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages"

# repeatMasker coordinates as GRangesList .RDS
rmsk_coords <- list.files(rmsk_coords_path, pattern = "rmsk\\.L1_and_LTRs\\.filtered.*\\.GRangesList\\.RDS", full.names = T)

######################################################## READ DATA
# read repeatMasker GRangesList 
rmsk_gr_list <- 
  map(rmsk_coords, readRDS) %>% 
  set_names(., c("whole", "removed_exons", "removed_genes"))

######################################################## MAIN CODE
# load bam file list in memory
bamfiles <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

### get count of reads, save summarizedExperiment as RDS
map(names(rmsk_gr_list), function(filter_name){
  
  # get features 
  rmsk_gr <- rmsk_gr_list[[filter_name]]
  
  # summarize overlaps - single end
  GenomicAlignments::summarizeOverlaps(features = rmsk_gr,
                                       reads = bamfiles,
                                       param = ScanBamParam(tagFilter = list("nM" = 0)),
                                       mode = "Union",
                                       singleEnd = TRUE,
                                       ignore.strand = TRUE) %>% 
    saveRDS(., file = file.path(outpath, str_c("rmsk.L1_and_LTRs", experiment_name, "20190306", filter_name, "SE.RDS", sep = ".")))
  
})
