### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
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
# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# set list of dataset names
dataset_name_list <- c("hamster_oocyte_Mov10l.deduplicated.smallRNAseq")

# get expression values for different datasets
purrr::map(dataset_name_list, function(dataset_name){
  
  ### set paths
  # mapped path
  base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", dataset_name)
  documentation_path <- file.path(base_path, "Data/Documentation")
  mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
  library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")
  
  # expression files path
  inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly"
  inpath <- file.path(inpath, "FLI_elements/LINE/expression", dataset_name, "expression_sum.RDS_files")
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
    dplyr::mutate(repName = "LINE1_full")
  
  # read library sizes
  library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 
  
  
  ### calculate RPMs for different subsets of repeats
  # create list of filtered table
  counts_tb_list <- list(LINE1_FLIs = expression_tb_full)
  
  # calculate RPMs in a loop
  RPM_summarized <- purrr::map(names(counts_tb_list), function(subset_name){
    
    # add category to expression table
    rpm_tb <- 
      counts_tb_list[[subset_name]] %>% 
      dplyr::group_by(sample_id, read_width, sense) %>% 
      dplyr::summarise(count = sum(count)) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::filter(!str_detect(sample_id, "HET")) %>% 
      dplyr::left_join(., library_size_tb, by = c("sample_id")) %>% 
      dplyr::mutate(library_size = library_size / 1e6, 
                    RPM = count / library_size) %>% 
      dplyr::mutate(sense = factor(sense, c("sense", "antisense"))) %>% 
      dplyr::arrange(sample_id, sense, read_width)
    
    # save as table
    rpm_tb %>% 
      tidyr::pivot_wider(., id_cols = "sample_id", names_from = c(read_width, sense), values_from = RPM, names_prefix = "r.", values_fill = 0) %T>% 
      readr::write_csv(file.path(outpath, str_c(subset_name, dataset_name, "RPM", "sense_antisense", "read_length.histogram", "20210517", "csv", sep = ".")))
    
    # return
    return(subset_name)
    
  })
  
  # return
  return(dataset_name)
  
})


