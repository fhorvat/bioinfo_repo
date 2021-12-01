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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/hamster_testis_Mov10l.8.5dpp.small_RNAseq"

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
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression"
gene_info_path <- file.path(gene_info_path, "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.geneInfo.csv")

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
  dplyr::bind_rows(.)

# # read gene info
# gene_info <- readr::read_csv(gene_info_path)
# 
# ### summarize gene info
# # get repName and repClass combinations
# gene_info_clean <- 
#   gene_info %>% 
#   dplyr::select(repName, repClass) %>% 
#   unique(.)

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
### all IAPs
# set name 
plot_name.1 <- "01.IAP_all"

# filter table
expression_tb.1 <- 
  expression_tb_full %>% 
  dplyr::filter(str_detect(repName, "^IAP"))


### LTR IAPs
# set name 
plot_name.2 <- "02.IAP_LTRs"

# filter table
expression_tb.2 <- 
  expression_tb_full %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter(!(str_detect(repName, "-int|_I$")))


### ints IAPs
# set name 
plot_name.3 <- "03.IAP_ints"

# filter table
expression_tb.3 <- 
  expression_tb_full %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter((str_detect(repName, "-int|_I$")))


### IAPLTR3/4, LTRs and ints
# set name 
plot_name.4 <- "04.IAPLTR3_4.LTRs_ints"

# filter table
expression_tb.4 <- 
  expression_tb_full %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter(repName %in% c("IAPLTR3", "IAPLTR3-int", "IAPLTR4", "IAPLTR4_I"))

# create list 
expression_list <- 
  list(expression_tb.1, expression_tb.2, expression_tb.3, expression_tb.4) %>% 
  set_names(., c(plot_name.1, plot_name.2, plot_name.3, plot_name.4))


### get expression 
# loop through expression
purrr::map(names(expression_list), function(subset_name){
  
  # get one expression table
  expression_tb <- expression_list[[subset_name]]
  
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
    dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_8.5dpp", "Mov10l_KO_8.5dpp")))
  
  # save as table
  count_tb %>% 
    tidyr::pivot_wider(., id_cols = "genotype_age", names_from = c(read_width, sense), values_from = RPM, names_prefix = "r.") %>% 
    dplyr::arrange(genotype_age) %T>%
    readr::write_csv(str_c(subset_name, "testis", "8.5dpp", "small_RNAseq", "RPM", "split_sense_antisense", "read_length.histogram.csv", sep = "."))
  
})

