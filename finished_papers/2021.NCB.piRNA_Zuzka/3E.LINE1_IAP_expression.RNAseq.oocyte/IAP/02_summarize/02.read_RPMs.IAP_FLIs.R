### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
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
# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# set list of dataset names
dataset_name_list <- c("hamster_oocyte_Mov10l.RNAseq", 
                       "hamster_testis_Mov10l.RNAseq", 
                       "hamster_testis_Mov10l.8.5dpp.RNAseq",
                       "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq")

# get expression values for different datasets
purrr::map(dataset_name_list, function(dataset_name){
  
  ### set paths
  # mapped path
  base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", dataset_name)
  documentation_path <- file.path(base_path, "Data/Documentation")
  sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
  mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
  library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")
  
  # expression files path
  inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly"
  inpath <- file.path(inpath, "FLI_elements/IAP/expression", dataset_name, "expression_sum.RDS_files")
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
    dplyr::bind_rows(.) %>% 
    dplyr::mutate(repName = "IAP_full")
  
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
  
  
  ### calculate RPMs for different subsets of repeats
  # create list of filtered table
  counts_tb_list <- list(IAP_FLIs = expression_tb_full)
  
  # calculate RPMs in a loop
  RPM_summarized <- purrr::map(names(counts_tb_list), function(subset_name){
    
    # add category to expression table
    rpm_tb <- 
      counts_tb_list[[subset_name]] %>% 
      dplyr::filter(!str_detect(sample_id, "HET")) %>% 
      dplyr::group_by(sample_id) %>% 
      dplyr::summarise(count = sum(count)) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::inner_join(., library_size_tb, by = c("sample_id")) %>% 
      dplyr::mutate(library_size = library_size / 1e6, 
                    RPM = count / library_size) %>% 
      dplyr::mutate(subset = subset_name) %>% 
      tidyr::pivot_wider(id_cols = subset, names_from = sample_id, values_from = RPM, names_prefix = "RPM.oocyte.")
    
    # save as table
    rpm_tb %T>% 
      readr::write_csv(file.path(outpath, str_c(subset_name, dataset_name, "RPM", "20201212", "csv", sep = ".")))
    
    # return
    return(subset_name)
    
  })
  
  # return
  return(dataset_name)
  
})


