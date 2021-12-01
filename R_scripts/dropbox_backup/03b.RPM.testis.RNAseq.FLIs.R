### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/library_composition")

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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/expression/hamster_testis_Mov10l.RNAseq"

# set outpath
outpath <- getwd()

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")

# gene info path
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE1"
gene_info_path <- file.path(gene_info_path, "LINE1.FLI_elements.csv")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

# calculated RPMs for other subsets
subset_rpms_path <- file.path(outpath, "LINE1s_subsets.testis.RNAseq.RPM.csv")

######################################################## READ DATA
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

# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- 
  readr::read_delim(library_size_path, delim = "\t") %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) 

# read calculated RPMs for other subsets
subset_rpms <- readr::read_csv(subset_rpms_path)

######################################################## MAIN CODE
### clean data
# filter sample table
sample_tb %<>% 
  dplyr::select(sample_id, genotype, age) %>% 
  tidyr::unite(genotype_age, genotype, age, sep = "_")

# summarize library size table
library_size_tb %<>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype_age) %>% 
  dplyr::summarise(library_size = sum(library_size)) 


### calculate RPMs for different subsets of repeats
# create list of filtered table
counts_tb_list <- list(LINE1_FLIs = expression_tb_full)

# calculate RPMs in loop
RPM_summarized <- purrr::map(names(counts_tb_list), function(subset_name){
  
  # add category to expression table
  rpm_tb <- 
    counts_tb_list[[subset_name]] %>% 
    dplyr::filter(!str_detect(sample_id, "HET")) %>% 
    dplyr::left_join(., sample_tb, by = "sample_id") %>% 
    dplyr::group_by(genotype_age) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., library_size_tb, by = c("genotype_age")) %>% 
    dplyr::mutate(library_size = library_size / 1e6, 
                  RPM = count / library_size) %>% 
    dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_13dpp", "Mov10l_KO_13dpp",
                                                                 "Mov10l_WT_21dpp", "Mov10l_KO_21dpp",
                                                                 "Mov10l_WT_adult", "Mov10l_KO_adult"))) %>%
    dplyr::arrange(genotype_age) %>% 
    dplyr::mutate(subset = subset_name) %>% 
    tidyr::pivot_wider(id_cols = subset, names_from = genotype_age, values_from = RPM, names_prefix = "RPM.testis.")
  
  # return
  return(rpm_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# join with already calculated RPMs for other subsets
RPM_summarized <- rbind(subset_rpms, RPM_summarized)

# save as table
readr::write_csv(RPM_summarized, file.path(outpath, str_c("LINE1s_subsets", "testis", "RNAseq", "RPM", "complete", "csv", sep = ".")))


