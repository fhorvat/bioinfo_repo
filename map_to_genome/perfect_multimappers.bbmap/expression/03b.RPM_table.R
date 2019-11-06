#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working dir
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
### IN AND OUT
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### COMMAND LINE ARGUMENTS
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)
dataset_path <- args$dataset_path
documentation_path <- args$documentation_path
feature_coordinates <- args$feature_coordinates

# reads stats path
read_stats_path <- file.path(dataset_path, "3_logs", "log.read_stats.txt")

# counts path
counts_path <- list.files(inpath, "LINE1_annotation.*\\.counts\\.txt$", full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read LINE1 table
line1_tb <-
  counts_tb %>%
  dplyr::select(rmsk_id = Geneid, seqnames = Chr, start = Start, end = End, strand = Strand, width = Length)

# read and clean stats
reads_stats <-
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>%
  mutate(library_size = round(library_size / 1E6, 6))

######################################################## MAIN CODE
# tidy table
rpm_tidy <-
  counts_tb %>%
  dplyr::select(-c(Chr:Length)) %>%
  dplyr::rename(rmsk_id = Geneid) %>%
  tidyr::gather(key = sample_id, value = counts, -rmsk_id) %>%
  mutate(sample_id = str_remove(sample_id, "\\.bam")) %>%
  dplyr::filter(sample_id != "s_GV_MILI_r1.PE") %>%
  left_join(., reads_stats, by = "sample_id") %>%
  mutate(rpm = round((counts / library_size), 4)) %>%
  dplyr::select(-c(counts, library_size))

# wide format
rpm_wide <-
  rpm_tidy %>%
  tidyr::spread(key = sample_id, value = rpm)

# average
rpm_average <-
  rpm_tidy %>%
  mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1|\\.[P,S]E|_old")) %>%
  group_by(rmsk_id, sample_id) %>%
  summarise(avg_rpm = mean(rpm) %>% round(., 3)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = sample_id, value = avg_rpm) %>%
  set_colnames(., c("rmsk_id", str_c(colnames(.)[-1], ".FPKM")))

# save
readr::write_csv(rpm_average, file.path(outpath, counts_path %>% basename(.) %>% str_replace(., "\\.counts\\.txt$", ".RPM.csv")))
