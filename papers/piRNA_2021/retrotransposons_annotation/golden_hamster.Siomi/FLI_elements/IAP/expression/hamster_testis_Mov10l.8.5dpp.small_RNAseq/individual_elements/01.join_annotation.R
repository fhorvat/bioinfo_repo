### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/expression/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq/individual_elements")

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

# set dataset name
dataset_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"

# annotation path
annotation_path <- "../../.."
annotation_path <- list.files(annotation_path, ".*\\.FLI_elements\\.csv", full.names = T)

# coverage path
coverage_path <- list.files(inpath, "\\.coverage_ratio\\.csv", full.names = T, recursive = T)

# expression path
expression_path <- list.files(inpath, "\\.RPM_expression\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read annotation
annotation_tb <- readr::read_csv(annotation_path)

# read coverage
coverage_tb <- readr::read_csv(coverage_path)

# read expression
expression_tb <- readr::read_csv(expression_path)

######################################################## MAIN CODE
# join all
annotation_tb %>% 
  dplyr::mutate(name = str_c(ltr_name, rmsk_id, sep = ".")) %>% 
  dplyr::left_join(., expression_tb, by = c("name" = "rmsk_id")) %>% 
  dplyr::left_join(., coverage_tb, by = "name") %>% 
  dplyr::select(-name) %T>%
  readr::write_csv(., file.path(outpath, basename(annotation_path) %>% str_replace(., "\\.csv$", str_c(".", dataset_name, ".expression_and_coverage.csv"))))
