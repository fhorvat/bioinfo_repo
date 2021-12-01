#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


## command line
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# experiment
experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)


### Base
# basepath
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression"


## Documentation
# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation"

# LINE1s coordinates
line1_path <- file.path(documentation_path, "2019_RNAi_piRNA_file3_L1_table.csv")


## Datasets
# datasets path
datasets_path <- file.path(base_path, "datasets", experiment, "Mapped/perfect_alignments.all_multimappers")

# get list of bam files
bam_path <- list.files(path = datasets_path, pattern = "*.bam$", full.names = T)

# summarized overlaps path
se_path <- list.files(inpath, ".*\\.se\\.RDS", full.names = T)

# reads stats path
read_stats_path <- file.path(datasets_path, "3_logs", "log.read_stats.txt")

# coverage table
coverage_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/element_coverage/2019_RNAi_piRNA_file3_L1_table.Dicer_Mili_KO.coverage.csv"

######################################################## READ DATA
# read LINE1 table
line1_tb <-
  readr::read_csv(line1_path) %>%
  dplyr::select_at(vars(-starts_with("X")))

# read and clean stats
reads_stats <- 
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  mutate(library_size = round(library_size / 1E6, 6))

# read coverage table
coverage_tb <- readr::read_csv(coverage_path)

# read summarizedOverlaps
counts_tb <- 
  readRDS(se_path) %>% 
  assay(.) %>% 
  as_tibble(., rownames = "rmsk_id")

######################################################## MAIN CODE
# tidy table
rpm_tidy <- 
  counts_tb %>% 
  tidyr::gather(key = sample_id, value = counts, -id) %>% 
  mutate(sample_id = str_remove(sample_id, "\\.bam")) %>% 
  dplyr::filter(sample_id != "s_GV_MILI_r1.PE") %>% 
  left_join(., reads_stats, by = "sample_id") %>%
  mutate(rpm = round((counts / library_size), 4)) %>% 
  dplyr::select(-c(counts, library_size))

# wide format
rpm_wide <- 
  rpm_tidy %>% 
  tidyr::spread(key = sample_id, value = rpm) %>% 
  right_join(., line1_coords, by = "id") %>% 
  dplyr::select_at(vars(id, seqnames:strand, repName, unique(rpm_tidy$sample_id))) 

# average
rpm_average <-
  rpm_tidy %>%
  mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1|\\.[P,S]E|GV_")) %>%
  group_by(id, sample_id) %>%
  summarise(avg_rpm = mean(rpm) %>% round(., 3)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = sample_id, value = avg_rpm) %>%
  left_join(., line1_coords, by = "id") %>%
  dplyr::select_at(vars(id, seqnames:strand, repName, starts_with("s_")))


### save
# open workbook and sheets
wb <- createWorkbook("RPM")
addWorksheet(wb, sheetName = "perfect_all_multimappers.RPM")
addWorksheet(wb, sheetName = "perfect_all_multimappers.mean_RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "perfect_reads.RPM", 
               x = rpm_wide %>% dplyr::select_at(vars(seqnames, start, end, strand, repName, starts_with("s_"))))

writeDataTable(wb = wb, 
               sheet = "perfect_reads.average_RPM", 
               x = rpm_average %>% dplyr::select_at(vars(seqnames, start, end, strand, repName, starts_with("s_"))))

# save workbook to disk
saveWorkbook(wb, file.path(results_path, str_c(basename(line1_path) %>% str_remove(., "\\.csv$"), experiment, "RPM.20191009.xlsx", sep = ".")), overwrite = TRUE)


### coverage table
# join with coverage table and save
line1_coverage_rpm <- 
  left_join(coverage_tb, rpm_wide %>% dplyr::select_at(vars("id", starts_with("s_"))), by = c("L1base_id" = "rmsk_id")) %>% 
  write_csv(., file.path(outpath, str_c(basename(results_path) %>% str_replace(., "\\.csv", ".RPM.csv"))))

