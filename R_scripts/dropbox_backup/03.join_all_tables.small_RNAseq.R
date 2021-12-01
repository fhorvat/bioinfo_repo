### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/expression/small_RNAseq")

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

# RPM histograms path
rpm_paths <- list.files(inpath, ".*read_length\\.histogram\\.[0-9]{8}\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read and join all tables
rpm_tb <- purrr::map(rpm_paths, function(path){
  
  # subset name
  subset_name <- 
    path %>% 
    basename(.) %>% 
    str_remove(., "\\.testis.*\\.csv") %>% 
    str_remove(., "^0[1-9]{1}\\.")
  
  # read table
  rpm <- 
    readr::read_csv(path) %>% 
    dplyr::mutate(subset = subset_name) %>% 
    dplyr::select(sample_id, subset, everything())
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# to long
rpm_tb_long <- 
  rpm_tb %>%
  dplyr::mutate(subset = factor(subset, levels = c("IAP_all", "IAP_LTRs", "IAP_ints", "IAPLTR3_4.LTRs_ints", "IAP_FLIs"))) %>% 
  tidyr::pivot_longer(cols = -c(sample_id, subset), names_to = "width", values_to = "RPM") %>% 
  tidyr::separate(col = width, into = c("width", "sense"), sep = "_") %>% 
  dplyr::mutate(width = width %>% str_remove(., "r\\.") %>% as.numeric(.))

### get histogram of read widths
# get only FLI elements
rpm_tb_hist <- 
  rpm_tb_long %>% 
  dplyr::filter(subset == "IAP_FLIs", 
                sense == "antisense") %>% 
  tidyr::pivot_wider(id_cols = sample_id, names_from = width, values_from = RPM, names_prefix = "antisense.r.") %>% 
  dplyr::mutate(age = str_extract(sample_id, "13dpp|21dpp|8.5_dpp|adult"), 
                genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_HET|Mov10l_KO"), 
                reseq = ifelse(str_detect(sample_id, "reseq"), T, F)) %>%
  dplyr::mutate(age = replace(age, age == "8.5_dpp", "8.5dpp"),
                age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp", "adult")),
                genotype = ifelse(reseq, str_c(genotype, ".reseq"), genotype),
                genotype = factor(genotype, levels = c("Mov10l_WT", "Mov10l_WT.reseq", "Mov10l_HET", "Mov10l_KO", "Mov10l_KO.reseq"))) %>% 
  dplyr::arrange(age, genotype)

# save for individual samples
rpm_tb_hist %>% 
  dplyr::select(sample_id, antisense.r.19:antisense.r.32) %T>%
  readr::write_csv(., file.path(outpath, str_c("IAPs_FLI.small_RNAseq", 
                                               str_c(19, "to", 32, "nt"), 
                                               "antisense", 
                                               "all_samples_RPM.20201102.csv", sep = ".")))

# calculate mean for genotype/age, save
rpm_tb_hist_mean <- 
  rpm_tb_long %>% 
  dplyr::filter(subset == "IAP_FLIs", 
                sense == "antisense") %>% 
  dplyr::mutate(age = str_extract(sample_id, "13dpp|21dpp|8.5_dpp|adult"), 
                genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_HET|Mov10l_KO"), 
                reseq = ifelse(str_detect(sample_id, "reseq"), T, F)) %>%
  dplyr::mutate(age = replace(age, age == "8.5_dpp", "8.5dpp"), 
                age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp", "adult")),
                genotype = ifelse(reseq, str_c(genotype, ".reseq"), genotype),
                genotype = factor(genotype, levels = c("Mov10l_WT", "Mov10l_WT.reseq", "Mov10l_HET", "Mov10l_KO", "Mov10l_KO.reseq"))) %>% 
  dplyr::arrange(age, genotype) %>% 
  dplyr::group_by(width, sense, age, genotype) %>% 
  dplyr::summarise(RPM = mean(RPM)) %>% 
  tidyr::pivot_wider(id_cols = c(age, genotype), names_from = width, values_from = RPM, names_prefix = "antisense.r.") %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_")

# save
rpm_tb_hist_mean %>% 
  readr::write_csv(., file.path(outpath, str_c("IAPs_FLI.small_RNAseq", 
                                               str_c(19, "to", 32, "nt"), 
                                               "antisense",
                                               "mean_RPM.20201102.csv", sep = ".")))

### filter by read length
# set read length
min_read <- 21
max_read <- 23

# long
rpm_filtered <- 
  rpm_tb_long %>% 
  dplyr::filter(width >= min_read, width <= max_read) %>% 
  dplyr::group_by(subset, sample_id, sense) %>% 
  dplyr::summarise(RPM = sum(RPM)) %>%
  dplyr::mutate(RPM = replace(RPM, is.na(RPM), 0)) %>% 
  dplyr::mutate(age = str_extract(sample_id, "GV|13dpp|21dpp|8.5_dpp|adult"), 
                genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_HET|Mov10l_KO"), 
                reseq = ifelse(str_detect(sample_id, "reseq"), T, F)) %>%
  dplyr::mutate(age = replace(age, age == "8.5_dpp", "8.5dpp"), 
                age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp", "adult")),
                genotype = ifelse(reseq, str_c(genotype, ".reseq"), genotype),
                genotype = factor(genotype, levels = c("Mov10l_WT", "Mov10l_WT.reseq", "Mov10l_HET", "Mov10l_KO", "Mov10l_KO.reseq")),
                sense = factor(sense, c("sense", "antisense"))) %>% 
  dplyr::arrange(sense, age, genotype)

# wide 
rpm_wide <- 
  rpm_filtered %>% 
  tidyr::pivot_wider(id_cols = subset, 
                     names_from = c(sample_id, sense), 
                     values_from = RPM, 
                     names_prefix = str_c("RPM.", min_read, "to", max_read, "nt.")) %T>%
  readr::write_csv(., file.path(outpath, str_c("IAPs_subsets.small_RNAseq", 
                                               str_c(min_read, "to", max_read, "nt"), 
                                               "all_samples_RPM.20201102.csv", 
                                               sep = ".")))

# mean by age/genotype
rpm_mean <- 
  rpm_filtered %>% 
  dplyr::group_by(subset, sense, age, genotype) %>% 
  dplyr::summarise(RPM = mean(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = subset,
                     names_from = c("age", "genotype", "sense"), 
                     values_from = "RPM", 
                     names_prefix = str_c("RPM.", min_read, "to", max_read, "nt.")) %T>%
  readr::write_csv(., file.path(outpath, str_c("IAPs_subsets.small_RNAseq", 
                                               str_c(min_read, "to", max_read, "nt"), 
                                               "mean_RPM.20201102.csv", 
                                               sep = ".")))

