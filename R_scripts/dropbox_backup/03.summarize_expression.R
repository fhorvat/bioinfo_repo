### INFO: 
### DATE: Sun Dec 13 03:23:07 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/expression/long_RNAseq")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

######################################################## READ DATA
# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

######################################################## MAIN CODE
# select relevant columns in classes table
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

# set list of dataset names
dataset_name_list <- c("hamster_oocyte_Mov10l.RNAseq", 
                       "hamster_testis_Mov10l.RNAseq", 
                       "hamster_testis_Mov10l.8.5dpp.RNAseq", 
                       "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq", 
                       "hamster_testis_Mov10l.0.5dpp.RNAseq", 
                       "hamster_MII_Piwil1.RNAseq")

# get expression values for different datasets
expression_tb_full <- purrr::map(dataset_name_list, function(dataset_name){
  
  ### set paths
  # mapped path
  base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", dataset_name)
  documentation_path <- file.path(base_path, "Data/Documentation")
  sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
  mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
  library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")
  
  # expression files path
  inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly"
  inpath <- file.path(inpath, "raw_rmsk.expression", "expression/long_RNAseq", dataset_name, "expression_sum.RDS_files")
  expression_path <- list.files(path = inpath, pattern = ".*\\.expression\\.RDS$", full.names = T)
  
  
  ### read data
  # read expression values
  expression_tb_full <- purrr::map(expression_path, function(path){
    
    # get sample name
    sample_name <- path %>% basename(.) %>% str_remove(., "\\.expression\\.RDS")
    
    # read table
    tmp_tb <- 
      path %>% 
      readRDS(.) %>% 
      dplyr::mutate(sample_id = sample_name)
    
    # return
    return(tmp_tb)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # read sample table
  sample_tb <- readr::read_csv(sample_table_path) 
  
  # read library sizes
  library_size_tb <- readr::read_delim(library_size_path, delim = "\t") 
  
  
  ### clean data
  # filter sample table
  sample_tb %<>% dplyr::select(sample_id, genotype) 
  
  # summarize library size table
  library_size_tb %<>% 
    dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
    dplyr::left_join(., sample_tb, by = "sample_id") %>% 
    dplyr::filter(!(sample_id %in% c("s_GV_Mov10l_KO_F66_r4.SE", "s_GV_Mov10l_HET_F53_r1.SE", "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE")))
  
  
  # add category to expression table, join with library size, get RPM
  rpm_tb <- 
    expression_tb_full %>% 
    dplyr::left_join(., classes_tb, by = "repName") %>% 
    dplyr::group_by(sample_id, category_I) %>% 
    # dplyr::group_by(sample_id, repName) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                  library_size = library_size / 1e6,
                  RPM = count / library_size) %>%
    dplyr::select(sample_id, category_I, RPM)
  
  # return
  return(rpm_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# get mean RPM per age/genotype/experiment
expression_tb_mean <- 
  expression_tb_full %>% 
  dplyr::mutate(tissue = str_extract(sample_id, "testis|GV|MII"), 
                tissue = replace(tissue, tissue == "GV", "oocyte"),
                tissue = replace(tissue, tissue == "MII", "egg"), 
                age = str_extract(sample_id, "GV|MII|0.5dpp|13dpp|21dpp|8.5dpp|adult|8.5"), 
                genotype = str_extract(sample_id, "Mov10l[1]*_WT|Mov10l[1]*_HET|Mov10l[1]*_KO|Piwil1_HET|Piwil1_KO") %>% 
                  str_replace(., "Mov10l_", "Mov10l1_")) %>%
  dplyr::mutate(age = replace(age, (sample_id %in% c("s_testis_Mov10l1_KO_8.5dpp_So820-M10_half_r1.SE", 
                                                     "s_testis_Mov10l1_KO_8.5dpp_So820-M10_r2.SE", 
                                                     "s_testis_Mov10l1_WT_8.5dpp_So802-M1_r1.SE", 
                                                     "s_testis_Mov10l1_WT_8.5dpp_So802-M3_r2.SE")), 
                              "8.5dpp.run_2")) %>% 
  dplyr::mutate(age = replace(age, age == "8.5", "8.5dpp")) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = c("oocyte", "egg", "testis")),
                age = factor(age, levels = c("GV", "MII", "0.5dpp", "8.5dpp", "8.5dpp.run_2", "13dpp", "21dpp", "adult")),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO", "Piwil1_HET", "Piwil1_KO"))) %>% 
  dplyr::arrange(tissue, age, genotype) %>% 
  dplyr::group_by(category_I, tissue, age, genotype) %>% 
  dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
  dplyr::ungroup(.) %>% 
  # dplyr::filter(!str_detect(genotype, "HET")) %>% 
  tidyr::pivot_wider(id_cols = c(category_I), names_from = c("tissue", "age", "genotype"), values_from = "RPM") %>% 
  dplyr::filter(!is.na(category_I)) %>% 
  dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))

### save 
# open workbook and sheets
wb <- createWorkbook("RPM")

# add sheet
addWorksheet(wb, sheetName = "RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "RPM", 
               x = expression_tb_mean)

# save workbook to disk
saveWorkbook(wb, file.path(outpath, str_c("all_LTR_classes.20210512", "long_RNAseq", "mean_RPM", "xlsx", sep = ".")), overwrite = TRUE)


