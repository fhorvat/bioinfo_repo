### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression")

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
line1_coords_path <- file.path(inpath, "Documentation", "L1s_nested_ours_20180516.ZJM.tidy.csv")

# datasets path
datasets_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/datasets" 

# results path
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/results"

# summarized overlaps path
se_path <- list.files(results_path, ".*se.RDS", full.names = T)

# reads stats path
read_stats_path <- list.files(datasets_path, pattern = "log\\..*\\.stats_and_tracks\\.csv", recursive = T, full.names = T)

# original table path
original_line1_path <- file.path(inpath, "Documentation", "L1s_nested_ours_20180517_PS_ZJM.xlsx")
  
######################################################## READ DATA
# read perfect stats
reads_stats <- 
  map(read_stats_path, function(path){
    
    # read and clean stats and tracks
    read_csv(path) %>%
      dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
      mutate(library_size = round(library_size / 1E6, 6), 
             experiment = str_remove(path, "\\/Mapped.*$") %>% basename(.) %>% str_remove(., "_[0-9]{4}.*")) %>%
             {if(str_detect(basename(path), "perfect")) mutate(., sample_id = str_c(sample_id, ".perfect")) else .}
    
  }) %>% 
  bind_rows(.)

# read Zoe's list of LINE1 full length elements
line1_coords <- read_csv(line1_coords_path)
original_line1_table <- openxlsx::read.xlsx(original_line1_path, sheet = "L1s_nested_ours_20180516_ZJM5")
  
# read summarizedOverlaps, join to one table
counts_tb <- map(se_path, function(path){
  
  readRDS(path) %>% 
    assay(.) %>% 
    as.tibble(., rownames = "gene_id")
  
}) %>% 
  purrr::reduce(., left_join, by = "gene_id")

######################################################## MAIN CODE
# tidy table
rpm_tidy <- 
  counts_tb %>% 
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>% 
  mutate(sample_id = str_remove(sample_id, "\\.bam")) %>% 
  left_join(., reads_stats, by = "sample_id") %>%
  mutate(rpm = round((counts / library_size), 4), 
         sample_id = str_c(sample_id, experiment, sep = ".")) %>% 
  dplyr::select(-c(counts, library_size, experiment))

# wide format
rpm_wide <- 
  rpm_tidy %>% 
  tidyr::spread(key = sample_id, value = rpm) %>% 
  right_join(., line1_coords, by = c("gene_id" = "id")) %>% 
  # dplyr::select_at(vars(seqnames:strand, repName, starts_with("s_"))) %>% 
  dplyr::select_at(vars(gene_id, contains("perfect")))

# # average
# rpm_average <- 
#   rpm_tidy %>% 
#   mutate(sample_id = str_remove_all(sample_id, "_r[0-9]{1,}|(?<=WT)[1-3]{1}|_BC[1-9]{1,}|(?<=WT)_1")) %>% 
#   group_by(gene_id, sample_id) %>% 
#   summarise(avg_rpm = mean(rpm) %>% round(., 3)) %>%
#   dplyr::ungroup(.) %>% 
#   tidyr::spread(key = sample_id, value = avg_rpm) %>% 
#   left_join(., line1_coords, by = c("gene_id" = "id")) %>% 
#   dplyr::select_at(vars(seqnames:strand, repName, starts_with("s_")))
  
### save
# open workbook and sheets
wb <- createWorkbook("RPM")
addWorksheet(wb, sheetName = "perfect_reads.RPM")
addWorksheet(wb, sheetName = "all_reads.RPM")
addWorksheet(wb, sheetName = "perfect_reads.average_RPM")
addWorksheet(wb, sheetName = "all_reads.average_RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "perfect_reads.RPM", 
               x = rpm_wide %>% 
                 dplyr::select_at(vars(seqnames, start, end, strand, repName, contains("perfect"))) %>% 
                 set_colnames(., str_remove(colnames(.), ".perfect")))

writeDataTable(wb = wb, 
               sheet = "all_reads.RPM", 
               x = rpm_wide %>% 
                 dplyr::select_at(vars(-contains("perfect"))))

writeDataTable(wb = wb, 
               sheet = "perfect_reads.average_RPM", 
               x = rpm_average %>% 
                 dplyr::select_at(vars(seqnames, start, end, strand, repName, contains("perfect"))) %>% 
                 set_colnames(., str_remove(colnames(.), ".perfect")))

writeDataTable(wb = wb, 
               sheet = "all_reads.average_RPM", 
               x = rpm_average %>% 
                 dplyr::select_at(vars(-contains("perfect"))))

# save workbook to disk
saveWorkbook(wb, file.path(results_path, "L1s_nested_ours_20180516.ZJM.WT_RNAseq.RPM.20190227.xlsx"), overwrite = TRUE)
