### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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

library(GEOquery)
library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment
experiment <- "Choi_2008_BiolReprod_GSE7774"

# base experiment path
base_experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays"

# experiment path 
experiment_path <- file.path(base_experiment_path, experiment)

# set working dir
setwd(experiment_path)

# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list .CEL files
cel_files <- list.files(inpath, pattern = "\\.CEL\\.gz", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
### experiment
# GEO accession
geo_accession <- str_remove(experiment, ".+_")

# get table from GEO, save
# geo_matrix <- getGEO(geo_accession, destdir = outpath)

anno_tb <- getGEO(file = "GPL1261.soft")

# tidy .cel files
cel_tb <- 
  tibble(cel_path = cel_files) %>% 
  dplyr::mutate(geo_accession = cel_path %>% basename(.) %>% str_remove(., "\\.CEL\\.gz"))

# tidy sample table
sample_table <- 
  geo_matrix %>% 
  .[[1]] %>% 
  pData(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(geo_accession, sample_desc_short = source_name_ch1, 
                strain = characteristics_ch1, sex = characteristics_ch1.1, 
                tissue = characteristics_ch1.2, age = characteristics_ch1.3) %>% 
  dplyr::mutate(stage = "newborn", 
                genotype = str_extract(sample_desc_short, "LHX8 Null|WT") %>% 
                  str_replace_all(., " ", "_")) %>% 
  dplyr::mutate_at(vars(strain, sex, tissue, age), list(~str_remove(., ".+: ") %>% tolower(.))) %>% 
  dplyr::mutate(tmp = str_c("s", stage, genotype, sep = "_")) %>% 
  dplyr::group_by(tmp) %>% 
  dplyr::mutate(sample_id = str_c(tmp, "_r_", 1:n())) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(sample_id, geo_accession, stage, genotype, tissue, strain, sex, age, sample_desc_short) %>% 
  dplyr::left_join(., cel_tb, by = "geo_accession") %T>% 
  readr::write_csv(., str_c(experiment, ".sampleTable.csv"))


### tidy annotation
# ## short annotation
# # get chip annotation from soft file
# chip_annotation <- getGEO(filename = "GPL1261.soft")
# 
# # get table, tidy
# chip_annotation_tb <-
#   eset_annotation@dataTable@table %>% 
#   as_tibble(.) %>% 
#   data.table::setnames(., colnames(.) %>% tolower(.) %>% str_replace_all(., " ", "_"))



