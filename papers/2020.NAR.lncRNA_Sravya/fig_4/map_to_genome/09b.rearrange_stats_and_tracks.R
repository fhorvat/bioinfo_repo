#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: 
### DATE: Thu Jul 11 17:26:04 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# stats and tracks path
tb_path <- list.files(inpath, "log.*.stats_and_tracks.csv")

######################################################## READ DATA
# read stats and tracks
stats_tb <- readr::read_csv(tb_path)

######################################################## MAIN CODE
# change columns
stats_tb %>% 
  dplyr::select(sample_id, raw.coverage, RPM_scaled.coverage, raw.individual_reads, 
                raw_input, mapped_total, unmapped_total,
                genome.mapped_minus_rDNA,
                genome.mapped_total, genome.uniquely_mapped, genome.mapped_to_multiple_loci, genome.unmapped_total, 
                total, rRNA, repeats, exon, other, 
                raw_input_rDNA = rDNA_45S.raw_input, total_mapped_rDNA = rDNA_45S.mapped_total, 
                uniquely_mapped_rDNA = rDNA_45S.uniquely_mapped, multi_mapped_rDNA = rDNA_45S.mapped_to_multiple_loci, unmapped_rDNA = rDNA_45S.unmapped) %>% 
  readr::write_csv(., file.path(outpath, basename(tb_path)))
