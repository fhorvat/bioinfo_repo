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
outpath <- getwd()

# find all tables with RPM values
rpm_tb_path <- list.files(inpath, ".*\\.RPM\\.complete\\.[0-6]{8}\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read and join table
rpm_tb <- purrr::map(rpm_tb_path, function(path){
  
  # read table
  rpm <- readr::read_csv(path)
  
  # change subset for small RNA-seq (wrong name)
  rpm %<>% 
    dplyr::mutate(subset = replace(subset, subset == "intact_guys", "IAP_FLIs"))
  
  # return
  return(rpm)
  
}) %>% 
  purrr::reduce(., left_join, by = "subset")

######################################################## MAIN CODE
# save table
readr::write_csv(rpm_tb, file.path(outpath, str_c("IAPs_subsets.long_RNAseq", 
                                                  "all_samples_RPM.20201102.csv", 
                                                  sep = ".")))

# calculate mean values
rpm_mean <- 
  rpm_tb %>% 
  dplyr::mutate(subset = factor(subset, levels = subset)) %>% 
  tidyr::pivot_longer(-subset, names_to = "sample_id", values_to = "RPM") %>% 
  dplyr::mutate(tissue = str_extract(sample_id, "testis|oocyte"), 
                age = str_extract(sample_id, "GV|13dpp|21dpp|8.5dpp|adult|8.5"), 
                genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_HET|Mov10l_KO")) %>%
  dplyr::mutate(age = replace(age, age == "8.5", "8.5dpp")) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = c("oocyte", "testis")),
                age = factor(age, levels = c("GV", "8.5dpp", "13dpp", "21dpp", "adult")),
                genotype = factor(genotype, levels = c("Mov10l_WT", "Mov10l_HET", "Mov10l_KO"))) %>% 
  dplyr::group_by(subset, tissue, age, genotype) %>% 
  dplyr::summarise(RPM = mean(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = subset, names_from = c("tissue", "age", "genotype"), values_from = "RPM")

# save
readr::write_csv(rpm_mean, file.path(outpath, str_c("IAPs_subsets.long_RNAseq", 
                                                    "mean_RPM.20201102.csv", 
                                                    sep = ".")))
