#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/counts/Dicer_Mili_KO/rmsk")

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

# experiment
experiment <- "Dicer_Mili_KO"

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/rmsk"

# LINE1s coordinates
line1_path <- file.path(documentation_path, "LINE1_annotation.rmsk.saf")

# basepath
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression"

# datasets path
datasets_path <- file.path(base_path, "datasets", experiment, "Mapped/perfect_alignments.all_multimappers")

# reads stats path
read_stats_path <- file.path(datasets_path, "3_logs", "log.read_stats.txt")

# counts path
counts_path <- list.files(inpath, "LINE1_annotation.*\\.counts\\.txt$", full.names = T)

# coverage path 
coverage_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/element_coverage.rmsk/Dicer_Mili_KO/LINE1.5k_to_7k.Dicer_Mili_KO.coverage.csv"

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

# read coverage table
coverage_tb <- readr::read_csv(coverage_path)

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
  mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1|\\.[P,S]E|GV_|_old")) %>%
  group_by(rmsk_id, sample_id) %>%
  summarise(avg_rpm = mean(rpm) %>% round(., 3)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = sample_id, value = avg_rpm) %>% 
  set_colnames(., c("rmsk_id", str_replace(colnames(.)[-1], "s_", "FPKM.")))

# join with coverage, save
coverage_rpm <- 
  coverage_tb %>% 
  dplyr::left_join(., rpm_average, by = "rmsk_id") %T>% 
  write.csv(., file.path(outpath, coverage_path %>% basename(.) %>% str_replace(., "coverage.csv", "coverage.RPM.csv")))
  


