### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/expression/long_RNAseq")

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
outpath <- file.path(getwd(), "joined_tables")

# find all tables with RPM values
rpm_tb_path <- list.files(inpath, "IAP_.*\\.RPM\\.[0-6]{8}\\.csv|IAPLTR3_4\\.LTRs_ints.*\\.csv", full.names = T, recursive = F)

######################################################## READ DATA
# read and join table
rpm_tb <- purrr::map(rpm_tb_path, function(path){
  
  # read table
  rpm <- 
    readr::read_csv(path) %>% 
    tidyr::pivot_longer(cols = -subset, names_to = "sample_id", values_to = "RPM") %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "RPM\\.oocyte\\.|RPM\\.testis\\."))
  
  # return
  return(rpm)
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# get long table
rpm_long <- 
  rpm_tb %>% 
  dplyr::mutate(subset = factor(subset, levels = c("IAP_all", "IAP_LTRs", "IAP_ints", "IAPLTR3_4.LTRs_ints", "IAP_FLIs"))) %>% 
  dplyr::mutate(tissue = str_extract(sample_id, "testis|GV"), 
                tissue = replace(tissue, tissue == "GV", "oocyte"), 
                age = str_extract(sample_id, "GV|13dpp|21dpp|8.5dpp|adult|8.5"), 
                genotype = str_extract(sample_id, "Mov10l[1]*_WT|Mov10l[1]*_HET|Mov10l[1]*_KO") %>% 
                  str_replace(., "Mov10l_", "Mov10l1_")) %>%
  dplyr::mutate(age = replace(age, (sample_id %in% c("s_testis_Mov10l1_KO_8.5dpp_So820-M10_half_r1.SE", 
                                                     "s_testis_Mov10l1_KO_8.5dpp_So820-M10_r2.SE", 
                                                     "s_testis_Mov10l1_WT_8.5dpp_So802-M1_r1.SE", 
                                                     "s_testis_Mov10l1_WT_8.5dpp_So802-M3_r2.SE")), 
                              "8.5dpp.run_2")) %>% 
  dplyr::mutate(age = replace(age, age == "8.5", "8.5dpp")) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = c("oocyte", "testis")),
                age = factor(age, levels = c("GV", "8.5dpp", "8.5dpp.run_2", "13dpp", "21dpp", "adult")),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"))) %>% 
  dplyr::arrange(subset, tissue, age, genotype)

# save wide table
rpm_wide <- 
  rpm_long %>% 
  tidyr::pivot_wider(id_cols = subset, values_from = "RPM", names_from = "sample_id") %T>% 
  readr::write_csv(., file.path(outpath, str_c("IAPs_subsets.long_RNAseq", 
                                               "all_samples_RPM.20201212.csv", 
                                               sep = ".")))

# calculate mean values
rpm_mean <- 
  rpm_long %>% 
  dplyr::group_by(subset, tissue, age, genotype) %>% 
  dplyr::summarise(RPM = mean(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = subset, names_from = c("tissue", "age", "genotype"), values_from = "RPM")

# save
readr::write_csv(rpm_mean, file.path(outpath, str_c("IAPs_subsets.long_RNAseq", 
                                                    "mean_RPM.20201212.csv", 
                                                    sep = ".")))
