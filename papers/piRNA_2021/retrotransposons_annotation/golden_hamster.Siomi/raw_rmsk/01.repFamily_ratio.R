### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/summary.repFamily")

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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/hamster_testis_Mov10l.small_RNAseq"

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# gene info path
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression"
gene_info_path <- file.path(gene_info_path, "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.geneInfo.csv")

# expression files path
expression_path <- file.path(inpath, "expression_sum.sense.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read expression values
expression_tb <- purrr::map(expression_path, function(path){
  
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

######################################################## MAIN CODE
### summarize gene info
# get repName and repClass combinations
gene_info_clean <-
  gene_info %>%
  dplyr::select(repName, repClass) %>%
  unique(.)

### filter library size
library_size_tb %<>% 
  dplyr::filter(str_detect(sample_id, "\\.24to31nt")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.24to31nt")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype, age) %>% 
  dplyr::summarise(library_size = sum(library_size))

### get expression 
# add category to expression table
count_tb <- 
  expression_tb %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::left_join(., gene_info_clean, by = "repName") %>% 
  dplyr::group_by(sample_id, repClass, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) 

# filter table
count_filt <- 
  count_tb %>% 
  dplyr::filter(read_width >= 24, read_width <= 31) %>% 
  dplyr::group_by(sample_id, repClass) %>%
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(sample_id, "HET")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(repClass, genotype, age) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = c("genotype", "age"), names_from = c(repClass), values_from = count) %>% 
  dplyr::left_join(., library_size_tb, by = c("genotype", "age"))


### calculate ratio of each class in whole library 
# retrotransposons
retrotrans_ratio <- 
  count_filt %>% 
  dplyr::mutate(LTR_all = ERV1 + ERVK + ERVL + LTR_other) %>% 
  dplyr::select(genotype, age, SINE, LINE, LTR_all, library_size) %>% 
  dplyr::mutate(SINE_percent = round((SINE / library_size), 3) * 100, 
                LINE_percent = round((LINE / library_size), 3) * 100, 
                LTR_percent = round((LTR_all / library_size), 3) * 100) %>% 
  dplyr::mutate(retro_percent = SINE_percent + LINE_percent + LTR_percent) %>% 
  dplyr::filter(genotype == "Mov10l_WT") %T>%
  readr::write_csv(., file.path(outpath, str_c("hamster_testis_Mov10l.small_RNAseq", "retrotrans", "24to31nt", "library_percent", "csv", sep = ".")))

# LTRs
ltrs_ratio <- 
  count_filt %>% 
  dplyr::mutate(LTR_all = ERV1 + ERVK + ERVL + LTR_other) %>% 
  dplyr::select(genotype, age, ERVK, ERVL, ERV1, LTR_other, LTR_all) %>% 
  dplyr::mutate(ERVK_percent = round((ERVK / LTR_all), 3) * 100,
                ERVL_percent = round((ERVL / LTR_all), 3) * 100, 
                ERV1_percent = round((ERV1 / LTR_all), 3) * 100,
                LTR_other_percent = round((LTR_other / LTR_all), 3) * 100) %>% 
  dplyr::filter(genotype == "Mov10l_WT") %T>%
  readr::write_csv(., file.path(outpath, str_c("hamster_testis_Mov10l.small_RNAseq", "LTRs", "24to31nt", "library_percent", "csv", sep = ".")))




