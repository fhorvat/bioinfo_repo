### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/expression/small_RNAseq")

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

# RPM histograms path
rpm_paths <- list.files(inpath, ".*read_length\\.histogram\\.[0-9]{8}\\.csv", full.names = T, recursive = T)

# samples list path
sample_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/datasets/Mov10l1_hamster.datasets.small_RNAseq.csv"

######################################################## READ DATA
# read and join all tables
rpm_tb <- purrr::map(rpm_paths, function(path){
  
  # subset name
  subset_name <- 
    path %>% 
    basename(.) %>% 
    str_remove(., "\\..*") 
  
  # read table
  rpm <- 
    readr::read_csv(path) %>% 
    dplyr::mutate(subset = subset_name) %>% 
    dplyr::select(sample_id, subset, everything())
  
}) %>% 
  dplyr::bind_rows(.)

# read samples table
sample_tb <- 
  readr::read_csv(sample_tb_path) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_HET|Mov10l_KO|Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = str_replace(genotype, "Mov10l_", "Mov10l1_"), 
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO")), 
                age = str_extract(sample_id, "13dpp|21dpp|8.5_dpp|adult|8.5dpp"), 
                age = replace(age, age == "8.5_dpp", "8.5dpp"),
                age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp", "adult"))) %>% 
  dplyr::filter((age == "21dpp" & experiment == "hamster_testis_Mov10l.small_RNAseq") | 
                  (age == "13dpp" & experiment == "hamster_testis_Mov10l.small_RNAseq.reseq") | 
                  (experiment == ("hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq")))

######################################################## MAIN CODE
# to long
rpm_tb_long <- 
  rpm_tb %>%
  dplyr::mutate(subset = factor(subset, levels = c("LINE1_all", "LINE1_Lx6_Lx5", "LINE1_FLIs"))) %>% 
  tidyr::pivot_longer(cols = -c(sample_id, subset), names_to = "width", values_to = "RPM") %>% 
  tidyr::separate(col = width, into = c("width", "sense"), sep = "_") %>% 
  dplyr::mutate(width = width %>% str_remove(., "r\\.") %>% as.numeric(.))


### get histogram of read widths
# apply to different subset
purrr::map(unique(rpm_tb_long$subset), function(element_subset){
  
  # get one element  
  rpm_tb_hist <- 
    rpm_tb_long %>% 
    dplyr::filter(subset == element_subset) %>% 
    dplyr::mutate(sense = factor(sense, levels = c("antisense", "sense"))) %>%
    dplyr::arrange(sense) %>% 
    tidyr::pivot_wider(id_cols = sample_id, names_from = c(width, sense), values_from = RPM, names_prefix = "r.") %>%
    dplyr::right_join(., sample_tb, by = "sample_id") %>% 
    dplyr::arrange(age, genotype)
  
  # save for individual samples
  rpm_tb_hist %>% 
    dplyr::select(sample_id, r.19_antisense:r.32_sense) %T>%
    readr::write_csv(x = ., file = file.path(outpath, str_c(element_subset, "small_RNAseq", 
                                                            str_c(19, "to", 32, "nt"), 
                                                            "sense_antisense", 
                                                            "all_samples_RPM.20210104.csv", sep = ".")))
  
  # calculate mean for genotype/age, save
  rpm_tb_hist_mean <- 
    rpm_tb_long %>% 
    dplyr::filter(subset == element_subset) %>% 
    dplyr::right_join(., sample_tb, by = "sample_id") %>% 
    dplyr::arrange(age, genotype) %>% 
    dplyr::group_by(width, sense, age, genotype) %>% 
    dplyr::summarise(RPM = mean(RPM)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(sense = factor(sense, levels = c("antisense", "sense"))) %>%
    dplyr::arrange(sense) %>% 
    tidyr::pivot_wider(id_cols = c(age, genotype), names_from = c(width, sense), values_from = RPM, names_prefix = "r.") %>% 
    tidyr::unite(col = "genotype_age", genotype, age, sep = "_")
  
  # save
  rpm_tb_hist_mean %>% 
    dplyr::select(genotype_age, r.19_antisense:r.32_sense) %T>%
    readr::write_csv(., file.path(outpath, str_c(element_subset, "small_RNAseq", 
                                                 str_c(19, "to", 32, "nt"), 
                                                 "sense_antisense",
                                                 "mean_RPM.20210104.csv", sep = ".")))
  
})


### summarize by read length
# set read length
read_length_list <- list(c(21, 23), c(24, 31), c(19, 32))

# loop through read lengths
purrr::map(read_length_list, function(read_length){
  
  # set min and max reads lengths
  min_read <- read_length[1]
  max_read <- read_length[2]
  
  # long
  rpm_filtered <- 
    rpm_tb_long %>% 
    dplyr::right_join(., sample_tb, by = "sample_id") %>% 
    dplyr::filter(width >= min_read, width <= max_read) %>% 
    dplyr::group_by(subset, sample_id, sense) %>% 
    dplyr::summarise(RPM = sum(RPM), 
                     age = unique(age), 
                     genotype = unique(genotype)) %>%
    dplyr::mutate(RPM = replace(RPM, is.na(RPM), 0)) %>% 
    dplyr::mutate(sense = factor(sense, c("sense", "antisense"))) %>% 
    dplyr::arrange(sense, age, genotype)
  
  # wide 
  rpm_wide <- 
    rpm_filtered %>% 
    tidyr::pivot_wider(id_cols = subset, 
                       names_from = c(sample_id, sense), 
                       values_from = RPM, 
                       names_prefix = str_c("RPM.", min_read, "to", max_read, "nt.")) %T>%
    readr::write_csv(., file.path(outpath, str_c("LINE1s_subsets.small_RNAseq", 
                                                 str_c(min_read, "to", max_read, "nt"), 
                                                 "all_samples_RPM.20210104.csv", 
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
    readr::write_csv(., file.path(outpath, str_c("LINE1s_subsets.small_RNAseq", 
                                                 str_c(min_read, "to", max_read, "nt"), 
                                                 "mean_RPM.20210104.csv", 
                                                 sep = ".")))
  
  # return 
  return(read_length)
  
})


