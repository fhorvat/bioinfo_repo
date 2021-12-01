### INFO: summarizes summarizedExperiments, calculates RPM
### DATE: Mon Feb 25 14:20:36 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
#wideScreen()

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

# clean rmsk path
rmsk_tidy_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/rmsk.L1_and_LTRs.filtered.20190306.csv"

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/developmental_stages.RNAseq.sample_table.20190225.csv"

# experiment
experiment_name <- "%EXPERIMENT"
experiment_name_short <- str_remove(experiment_name, "_.*")

# perfect read stats path
read_stats_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment_name, "Data/Mapped/STAR_mm10/4_count_perfect_reads")

# read stats list 
read_stats_list <- list.files(read_stats_path, pattern = ".*\\.read_stats.txt", full.names = T)

# summarizedExperiments list
se_path_list <- list.files(inpath, str_c("rmsk\\.L1_and_LTRs\\.", experiment_name, ".*\\.SE\\.RDS"))

######################################################## READ DATA
# read shorter list of retrotransposons from repeatMasker
rmsk_tidy <- 
  read_csv(rmsk_tidy_path) %>% 
  as.data.table(.)

# read sample table
sample_table <- read_csv(sample_table_path)

# read and sum read stats
read_stats <- 
  map(read_stats_list, function(path){
    
    # read and clean
    read_delim(path, delim = "\t") %>%
      dplyr::mutate(genome.mapped_minus_rDNA = ifelse(str_detect(sample_id, "\\.PE"), 
                                                      round(genome.mapped_minus_rDNA / 2), 
                                                      genome.mapped_minus_rDNA)) %>%
      dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
      mutate(library_size = round(library_size / 1E6, 6))
    
  }) %>% 
  bind_rows(.) %>% 
  as.data.table(.)

# read 3 different summarizedExperiments 
se_list <- 
  map(se_path_list, readRDS) %>% 
  set_names(., str_extract(se_path_list, "whole|removed_exons|removed_genes"))

######################################################## MAIN CODE
# filter sample table
sample_table_filt <- 
  sample_table %>% 
  dplyr::filter(experiment == experiment_name) %>%
  dplyr::select("sample_id", "stage") %>% 
  as.data.table(.)

### get average RPM for 3 different filtered repeatMaskers (whole, filtered exons, filtered genes)
# loop through 3 summarizedExperiments
rmsk_avgRPM_list <- map(names(se_list), function(se_name){
  
  # get counts
  avg_rpm_tb <- 
    se_list[[se_name]] %>% 
    assay(.) %>% 
    set_colnames(., str_remove(colnames(.), "\\.genome.*")) %>% 
    data.table::as.data.table(., keep.rownames = "rmsk_id") %>% 
    melt(., id.vars = "rmsk_id", variable.name = "sample_id",  value.name = "counts") %>% 
    .[read_stats, on = "sample_id", `:=`(rpm = round((counts / i.library_size), 3))] %>% 
    .[sample_table_filt, on = "sample_id", `:=`(stage = i.stage)] %>% 
    .[, list(avg_rpm = mean(rpm)), by = .(rmsk_id, stage)] %>% 
    .[rmsk_tidy[, c("rmsk_id", "repName", "repClass", "repFamily")], on = "rmsk_id", 
      `:=`(repName = i.repName, repClass = i.repClass, repFamily = i.repFamily)] %>% 
    .[]
  # dcast(., gene_id ~ stage, value.var = "avg_rpm")
  
  ### sum for different elements 
  # LINE 1
  line1_all <- 
    avg_rpm_tb[repFamily == "L1", 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repClass, stage)] %>% 
    .[, list(class = "LINE1_all", stage, sum_avg_rpm)]
  
  # Line1 subfamilies = L1Md_T L1Md_A L1Md_Gf L1Md_F2 L1_Mus2 L1_Mus1 L1Md_F
  line1_subfamilies <- 
    avg_rpm_tb[repName %in% c("L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repName, stage)]  %>% 
    .[, list(class = repName, stage, sum_avg_rpm)]
  
  # IAP all
  iap_all <- 
    avg_rpm_tb[str_detect(repName, "IAP"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repFamily, stage)] %>% 
    .[, list(class = "IAP_all", stage, sum_avg_rpm)]
  
  # IAPEz-Int with LTRs IAPLTR1_Mm
  iap_int_ltr <- 
    avg_rpm_tb[repName %in% c("IAPEy-int", "IAPLTR1_Mm"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repFamily, stage)] %>% 
    .[, list(class = "IAPEy_int_IAPLTR1_Mm", stage, sum_avg_rpm)]
  
  # MuERV-int with MT2 LTRs = MERVL-int + MT2_Mm
  mervl_mt2mm <- 
    avg_rpm_tb[repName %in% c("MERVL-int", "MT2_Mm"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repFamily, stage)] %>% 
    .[, list(class = "MERVL_int_MT2_Mm", stage, sum_avg_rpm)]
  
  # MTA all LTRs and ints
  mta_int_ltr <- 
    avg_rpm_tb[str_detect(repName, "MTA"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repFamily, stage)] %>% 
    .[, list(class = "MTA_int_LTR", stage, sum_avg_rpm)]
  
  # ORR1A all ints and LTR subgroups
  orr1a_int_ltr <- 
    avg_rpm_tb[str_detect(repName, "ORR1A"), 
               list(sum_avg_rpm = sum(avg_rpm)), 
               by = .(repFamily, stage)] %>% 
    .[, list(class = "ORR1A_int_LTR", stage, sum_avg_rpm)]
  
  ### bind all and save
  # bind rows to one table
  retro_all_tb <- 
    rbind(line1_all, line1_subfamilies, 
          iap_all, iap_int_ltr, 
          mervl_mt2mm, mta_int_ltr, 
          orr1a_int_ltr) %>% 
    .[, experiment := experiment_name_short] %>% 
    .[]
  
}) %>% 
  set_names(., names(se_list)) %T>% 
  saveRDS(., file = file.path(outpath, str_c("rmsk.L1_and_LTRs", experiment_name, "20190306", "avgRPM.RDS", sep = ".")))

