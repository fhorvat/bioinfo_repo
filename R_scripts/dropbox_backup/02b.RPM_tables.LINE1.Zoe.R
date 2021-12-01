#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
# set experiment
experiment <- "Dicer_SOM"
single_end <- TRUE

# set working dir
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/FPKM_tables/", experiment))

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

# list of LINE1 full length elements path
line1_coords_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

# datasets path
datasets_path <- file.path(inpath, "../../datasets", experiment, "Mapped/mm10_masked")

# results path
results_path <- outpath

# summarized overlaps path
se_path <- list.files(results_path, "L1s_nested_ours.*\\.se\\.RDS", full.names = T)

# reads stats path
read_stats_path <- list.files(datasets_path, pattern = "log.read_stats.txt", recursive = T, full.names = T)

# original table path
original_line1_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180517_PS_ZJM.xlsx"

######################################################## READ DATA
# read and clean stats and tracks
reads_stats <- 
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  mutate(library_size = round(library_size / 1E6, 6))

# read tidy Zoe's list of LINE1 full length elements
line1_coords <- read_csv(line1_coords_path)

# read original Zoe's table
line1_original_table <- 
  openxlsx::read.xlsx(original_line1_path, sheet = 2) %>% 
  set_colnames(., make.unique(colnames(.))) %>% 
  as_tibble(.)

# read summarizedOverlaps
counts_tb <- 
  readRDS(se_path) %>% 
    assay(.) %>% 
    as_tibble(., rownames = "id")

######################################################## MAIN CODE
# tidy table
rpm_tidy <- 
  counts_tb %>% 
  tidyr::gather(key = sample_id, value = counts, -id) %>% 
  mutate(sample_id = str_remove(sample_id, "\\.perfect\\.bam")) %>% 
  dplyr::filter(sample_id != "s_MII_B6_WT_1.SE") %>% 
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
addWorksheet(wb, sheetName = "perfect_reads.RPM")
addWorksheet(wb, sheetName = "perfect_reads.average_RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "perfect_reads.RPM", 
               x = rpm_wide %>% dplyr::select_at(vars(seqnames, start, end, strand, repName, starts_with("s_"))))

writeDataTable(wb = wb, 
               sheet = "perfect_reads.average_RPM", 
               x = rpm_average %>% dplyr::select_at(vars(seqnames, start, end, strand, repName, starts_with("s_"))))

# save workbook to disk
saveWorkbook(wb, file.path(results_path, str_c("L1s_nested_ours_20180516.ZJM.tidy.", experiment, ".RPM.20190929.xlsx")), overwrite = TRUE)


### original table
# join with original table and save
line1_original_table_rpm <- 
  left_join(line1_original_table, 
            rpm_wide %>% dplyr::select_at(vars("id", starts_with("s_"))), 
            by = c("uniqName" = "id")) %>% 
  write_csv(., file.path(results_path, str_c("L1s_nested_ours_20180517_PS_ZJM.original.", experiment, ".RPM.csv")))

# join with original table average RPM and save
line1_original_table_average_rpm <- 
  left_join(line1_original_table, 
            rpm_average %>% dplyr::select_at(vars("id", starts_with("s_"))), 
            by = c("uniqName" = "id")) %>% 
  write_csv(., file.path(results_path, str_c("L1s_nested_ours_20180517_PS_ZJM.original.", experiment, ".average_RPM.csv")))
