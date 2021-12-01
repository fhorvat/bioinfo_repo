### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/library_composition/antisense_reads_width_histogram.8.5dpp")

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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/expression/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "library_sizes.txt")

# gene info path
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP"
gene_info_path <- file.path(gene_info_path, "IAP.FLI_elements.csv")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

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
  dplyr::mutate(repName = "IAP_full")

# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 

# filter and summarize library size
library_size_tb %<>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") %>% 
  dplyr::group_by(genotype_age) %>% 
  dplyr::summarise(library_size = sum(library_size))

# filter sample table
sample_tb %<>% 
  dplyr::select(sample_id, genotype, age) %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") 

######################################################## MAIN CODE
### full length guys
# set name
plot_name <- "05.intact_guys"

# filter table
expression_tb <- expression_tb_full


### get expression 
# add category to expression table
count_tb <- 
  expression_tb %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::group_by(sample_id, read_width, sense) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(sample_id, "HET")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype_age, sense, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_size_tb, by = c("genotype_age")) %>% 
  dplyr::mutate(library_size = library_size / 1e6, 
                RPM = count / library_size) %>% 
  dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_8.5", "Mov10l_KO_8.5")))

# save as table
count_tb %>% 
  tidyr::pivot_wider(., id_cols = "genotype_age", names_from = c(read_width, sense), values_from = RPM, names_prefix = "r.") %>% 
  dplyr::arrange(genotype_age) %T>%
  readr::write_csv(str_c(plot_name, "testis", "8.5dpp", "small_RNAseq", "RPM", "split_sense_antisense", "read_length.histogram.csv", sep = "."))
