### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/whole_rmsk/expression/hamster_testis_Mov10l.small_RNAseq")

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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 

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
### filter library size
library_size_tb %<>% 
  dplyr::filter(str_detect(sample_id, "\\.24to31nt")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.24to31nt")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype, age) %>% 
  dplyr::summarise(library_size = sum(library_size))

### get expression 
# add category to expression table, join with library size, get RPM, join with sample table, summarize per genotype
count_tb <- 
  expression_tb %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::group_by(sample_id, rmsk_name, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) 
  
# filter table
count_filt <- 
  count_tb %>% 
  dplyr::filter(read_width >= 24, read_width <= 31) %>% 
  dplyr::group_by(sample_id, rmsk_name) %>%
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(sample_id, "HET")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(rmsk_name, genotype, age) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = c("genotype", "age"), names_from = c(rmsk_name), values_from = count) %>% 
  dplyr::left_join(., library_size_tb, by = c("genotype", "age"))

### calculate ratio of each class in whole library 
# retrotransposons
retrotrans_ratio <- 
  count_filt %>% 
  dplyr::mutate(LTR = ERV1 + ERVK + ERVL + LTR_other) %>% 
  dplyr::select(genotype, age, SINE, LINE, LTR, library_size) %>% 
  dplyr::mutate(SINE = round((SINE / library_size), 3) * 100, 
                LINE = round((LINE / library_size), 3) * 100, 
                LTR = round((LTR / library_size), 3) * 100) %>% 
  dplyr::mutate(retro = SINE + LINE + LTR) %>% 
  dplyr::filter(genotype == "Mov10l_WT")

# LTRs
ltrs_ratio <- 
  count_filt %>% 
  dplyr::mutate(LTR = ERV1 + ERVK + ERVL + LTR_other) %>% 
  dplyr::select(genotype, age, ERVK, ERVL, ERV1, LTR_other, LTR) %>% 
  dplyr::mutate(ERVK = round((ERVK / LTR), 3) * 100,
                ERVL = round((ERVL / LTR), 3) * 100, 
                ERV1 = round((ERV1 / LTR), 3) * 100,
                LTR_other = round((LTR_other / LTR), 3) * 100) %>% 
  dplyr::filter(genotype == "Mov10l_WT")


# open workbook and sheets
wb <- createWorkbook("RPM")

# add sheet
addWorksheet(wb, sheetName = "RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "RPM", 
               x = rpm_filt)

# save workbook to disk
saveWorkbook(wb, file.path(outpath, str_c("rmsk_repFamily.200730", "small_RNA.testis", "24to31nt", "mean_RPM", "xlsx", sep = ".")), overwrite = TRUE)


