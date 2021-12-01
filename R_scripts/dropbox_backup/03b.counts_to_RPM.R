#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts to RPM table 
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
feature_coordinates <- args$feature_coordinates

experiment='Manakov_2015_CellRep_GSE70731'
single_end=TRUE
threads=1
dataset_path='/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Manakov_2015_CellRep_GSE70731/Mapped/perfect_alignments.all_multimappers/4_merged_replicates'
feature_coordinates='/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/rmsk/LINE1_annotation.rmsk.GRanges.RDS'



### OTHER 
# summarizedOverlaps path
se_path <- list.files(inpath, pattern = "LINE1_annotation\\..*\\.se\\.RDS")

# coverage table
coverage_path <- list.files(inpath, pattern = "LINE1_annotation\\..*\\.coverage\\.csv")
  
# reads stats path
read_stats_path <- file.path(dataset_path, "../3_logs/log.read_stats.txt")


######################################################## READ DATA
# read LINE1 coordinates
line1_gr <- readRDS(feature_coordinates)

# read and clean stats
reads_stats <- 
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1|\\.[P,S]E|GV_")) %>%
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(library_size = sum(library_size), 
                   group_size = n()) %>% 
  dplyr::ungroup(.) %>% 
  mutate(library_size = round(library_size / 1E6, 6))

# read coverage table
coverage_tb <- 
  readr::read_csv(coverage_path) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

# read summarizedOverlaps
counts_tb <- 
  readRDS(se_path) %>% 
  assay(.) %>% 
  as_tibble(., rownames = "rmsk_id")

######################################################## MAIN CODE
# tidy table
rpm_tidy <- 
  counts_tb %>% 
  tidyr::gather(key = sample_id, value = counts, -rmsk_id) %>% 
  mutate(sample_id = str_remove(sample_id, "\\.bam")) %>% 
  left_join(., reads_stats, by = "sample_id") %>%
  mutate(rpm = round((counts / library_size), 4), 
         rpm = round((rpm / group_size), 3)) %>% 
  dplyr::select(-c(counts, library_size, group_size))

# wide format 
rpm_wide <- 
  rpm_tidy %>% 
  tidyr::spread(key = sample_id, value = rpm) %>% 
  set_colnames(., c("rmsk_id", str_c("FPKM", colnames(.)[-1], sep = ".")))

# join with coverage table and save
line1_coverage_rpm <- 
  coverage_tb %>% 
  left_join(., rpm_wide, by = "rmsk_id") %>% 
  write_csv(., file.path(outpath, str_c((feature_coordinates %>% basename(.) %>% str_remove(., "\\.GRanges\\.RDS$")), 
                                        experiment, 
                                        "coverage.FPKM.csv", 
                                        sep = ".")))

