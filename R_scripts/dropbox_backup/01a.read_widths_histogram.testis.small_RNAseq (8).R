### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/expression/small_RNAseq/hamster_testis_Mov10l.small_RNAseq.reseq")

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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly"
inpath <- file.path(inpath, "raw_rmsk.expression/hamster_testis_Mov10l.small_RNAseq.reseq")

# set outpath
outpath <- getwd()

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

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

# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 

######################################################## MAIN CODE
### prepare data
# filter and tidy sample table
sample_tb %<>%
  dplyr::select(sample_id, genotype, age)

# filter and tidy library size table
library_size_tb %<>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_")


### calculate RPMs for different subsets of repeats
# create list of filtered table
counts_tb_list <- list(IAP_all = expression_tb_full %>% 
                         dplyr::filter(str_detect(repName, "^IAP")), 
                       IAP_LTRs = expression_tb_full %>% 
                         dplyr::filter(str_detect(repName, "^IAP")) %>% 
                         dplyr::filter(!(str_detect(repName, "-int|_I$"))), 
                       IAP_ints = expression_tb_full %>% 
                         dplyr::filter(str_detect(repName, "^IAP")) %>% 
                         dplyr::filter((str_detect(repName, "-int|_I$"))), 
                       IAPLTR3_4.LTRs_ints = expression_tb_full %>% 
                         dplyr::filter(str_detect(repName, "^IAP")) %>% 
                         dplyr::filter(repName %in% c("IAPLTR3", "IAPLTR3-int", "IAPLTR4", "IAPLTR4_I")))

# calculate RPMs in a loop
RPM_summarized <- purrr::map(names(counts_tb_list), function(subset_name){
  
  # add category to expression table
  rpm_tb <- 
    counts_tb_list[[subset_name]] %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
    dplyr::group_by(sample_id, read_width, sense) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::filter(!str_detect(sample_id, "HET")) %>% 
    dplyr::left_join(., library_size_tb, by = c("sample_id")) %>% 
    dplyr::mutate(library_size = library_size / 1e6, 
                  RPM = count / library_size) %>% 
    dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_13dpp", "Mov10l_KO_13dpp", "Mov10l_KO_21dpp")), 
                  sense = factor(sense, c("sense", "antisense"))) %>% 
    dplyr::arrange(sample_id, sense, read_width)
  
  # save as table
  rpm_tb %>% 
    tidyr::pivot_wider(., id_cols = "sample_id", names_from = c(read_width, sense), values_from = RPM, names_prefix = "r.", values_fill = 0) %T>% 
    readr::write_csv(str_c(subset_name, "testis", "small_RNAseq.reseq", "RPM", "sense_antisense", "read_length.histogram", "20201102", "csv", sep = "."))
  
})
  
