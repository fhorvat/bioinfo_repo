#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
# set experiment
experiment <- "hamster_oocyte_Mov10l.RNAseq"
single_end <- TRUE

# set working dir
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/FPKM_tables/L1_joined_rmskID", experiment))

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

# documentation path 
documentation_path <- file.path(inpath, "../../../Documentation/L1_joined_rmskID/mesAur1.L1_joined_rmskID_20190929.FH")

# list of LINE1 full length elements path
line1_coords_path <- file.path(documentation_path, "LINE1_whole.joined_rmsk_ID.csv")

# datasets path
datasets_path <- file.path(inpath, "../../../datasets", experiment, "Mapped/mesAur1.L1_joined_rmskID_20190929.FH")

# results path
results_path <- outpath

# summarized overlaps path
se_path <- file.path(results_path, str_c("LINE1_whole.joined_rmsk_ID.", experiment, ".se.RDS"))

# reads stats path
read_stats_path <- list.files(datasets_path, pattern = "log.read_stats.txt", recursive = T, full.names = T)

######################################################## READ DATA
# read and clean stats and tracks
reads_stats <- 
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  mutate(library_size = round(library_size / 1E6, 6))

# read tidy list of LINE1 full length elements
line1_coords <- 
  read_csv(line1_coords_path) %>% 
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
  mutate(sample_id = str_remove(sample_id, "\\.perfect\\.bam")) %>% 
  left_join(., reads_stats, by = "sample_id") %>%
  mutate(rpm = round((counts / library_size), 4)) %>% 
  dplyr::select(-c(counts, library_size))

# wide format
rpm_wide <- 
  rpm_tidy %>% 
  tidyr::spread(key = sample_id, value = rpm) %>% 
  right_join(., line1_coords, by = "rmsk_id") %>% 
  dplyr::select_at(vars(rmsk_id, seqnames:strand, repName, overlaps_ensembl_id, unique(rpm_tidy$sample_id))) 

# average
rpm_average <-
  rpm_tidy %>%
  mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1|\\.[P,S]E|GV_")) %>%
  group_by(rmsk_id, sample_id) %>%
  summarise(avg_rpm = mean(rpm) %>% round(., 3)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = sample_id, value = avg_rpm) %>%
  left_join(., line1_coords, by = "rmsk_id") %>%
  dplyr::select_at(vars(rmsk_id, seqnames:strand, repName, overlaps_ensembl_id, starts_with("s_")))

### save
# open workbook and sheets
wb <- createWorkbook("RPM")
addWorksheet(wb, sheetName = "perfect_reads.RPM")
addWorksheet(wb, sheetName = "perfect_reads.average_RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "perfect_reads.RPM", 
               x = rpm_wide %>% dplyr::select_at(vars(rmsk_id, seqnames, start, end, strand, repName, overlaps_ensembl_id, starts_with("s_"))))

writeDataTable(wb = wb, 
               sheet = "perfect_reads.average_RPM", 
               x = rpm_average %>% dplyr::select_at(vars(rmsk_id, seqnames, start, end, strand, repName, overlaps_ensembl_id, starts_with("s_"))))

# save workbook to disk
saveWorkbook(wb, file.path(results_path, se_path %>% basename(.) %>% str_replace(., ".se.RDS", ".RPM.xlsx")), overwrite = TRUE)
