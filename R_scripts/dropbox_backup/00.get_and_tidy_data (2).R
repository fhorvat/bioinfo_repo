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
experiment <- "Joshi_2007_BMCDevBiol_GSE5558"

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

######################################################## READ DATA

######################################################## MAIN CODE
### experiment
# GEO accession
geo_accession <- str_remove(experiment, ".+_")

# get table from GEO, save
geo_matrix <- getGEO(geo_accession, destdir = outpath)

# tidy sample table
sample_table <- 
  geo_matrix %>% 
  .[[1]] %>% 
  pData(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(geo_accession, sample_desc_short = source_name_ch1, sample_desc_long = characteristics_ch1) %>% 
  dplyr::mutate(stage = str_extract(sample_desc_short, "E17.5|E14.5|E12.5|newborn"), 
                genotype = str_extract(sample_desc_short, "FIGLA|Wild-type") %>% 
                  str_replace_all(., c("FIGLA" = "Figla_Null", "Wild-type" = "WT")),
                tissue = "ovary", 
                replicate = str_extract(sample_desc_short, "replicate [0-9]+$") %>% 
                  str_replace(., "replicate ", "r_"), 
                sample_id = str_c("s", stage, genotype, replicate, sep = "_")) %>% 
  dplyr::select(sample_id, geo_accession, stage, genotype, tissue) %T>% 
  readr::write_csv(., str_c(experiment, ".sampleTable.csv"))


### tidy annotation
# ## short annotation
# # get chip annotation from soft file
# chip_annotation <- getGEO(filename = "GPL3771.soft")
# 
# # get table, tidy
# chip_annotation_tb <-
#   eset_annotation@dataTable@table %>% 
#   as_tibble(.) %>% 
#   data.table::setnames(., colnames(.) %>% tolower(.) %>% str_replace_all(., " ", "_"))



