### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV.repeatMasker/annotation")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# annotation table path
annotation_path <- file.path(inpath, "MYSERV6-int.longer_than_4kb.csv")

######################################################## READ DATA
# read annotation
annotation_tb <- readr::read_csv(annotation_path)

######################################################## MAIN CODE
# set list of dataset names
dataset_name_list <- c(
  # "hamster_oocyte_Mov10l.RNAseq", 
  # "hamster_testis_Mov10l.RNAseq", 
  # "hamster_testis_Mov10l.8.5dpp.RNAseq",
  "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq")

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
  inpath <- file.path(inpath, "MYSERV6", "expression/long_RNAseq", dataset_name, "expression_sum.RDS_files")
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
  
  ### calculate RPMs
  # calculate RPM
  rpm_tb <- 
    expression_tb_full %>% 
    dplyr::filter(!str_detect(sample_id, "HET")) %>% 
    dplyr::inner_join(., library_size_tb, by = c("sample_id")) %>% 
    dplyr::mutate(library_size = library_size / 1e6, 
                  RPM = count / library_size) %>% 
    dplyr::select(rmsk_coords, rmsk_id, sample_id, RPM)

  # return
  return(rpm_tb)
  
}) %>% 
  dplyr::bind_rows(.) 

# get mean RPM per age/genotype/experiment
expression_tb_mean <- 
  expression_tb_full %>% 
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
  dplyr::arrange(tissue, age, genotype) %>% 
  dplyr::group_by(rmsk_coords, tissue, age, genotype) %>% 
  dplyr::summarise(RPM = mean(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = c(rmsk_coords), names_from = c("tissue", "age", "genotype"), values_from = "RPM")

# add to annotation table
annotation_rpm <- 
  annotation_tb %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::left_join(., expression_tb_mean, by = c("coordinates" = "rmsk_coords")) %>% 
  dplyr::mutate_at(vars(contains("Mov10l1")), ~replace(., is.na(.), 0)) %T>%
  readr::write_csv(., file = file.path(outpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.csv"))


