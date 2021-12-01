### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE/brain/bigWig")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# read metadata table
table_path <- list.files(path = outpath, pattern = "metadata.tsv", full.names = T)
metadata <- readr::read_delim(file = table_path, delim = "\t")

######################################################## MAIN CODE
# create rename file
rename_out <- 
  metadata %>% 
  dplyr::select(file_accession = `File accession`, format = `File format`, 
                experiment_accesion = `Experiment accession`, sample_derived = `Derived from`, 
                output_type = `Output type`, biosample = `Biosample term name`, 
                assembly = Assembly) %>% 
  dplyr::mutate(strand = ifelse(str_detect(output_type, "plus"), "plus", "minus"), 
                read_group = ifelse(str_detect(output_type, "unique"), "unique", "all"), 
                experiment_group = data.table::rleid(experiment_accesion)) %>% 
  dplyr::arrange(experiment_accesion, sample_derived) %>% 
  base::split(., .$experiment_group) %>% 
  lapply(X = ., FUN = function(X) X %>% mutate(sample_group = data.table::rleid(sample_derived))) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::arrange(experiment_accesion, sample_derived, read_group, strand) %>% 
  dplyr::select(-c(output_type, experiment_accesion, sample_derived)) %>% 
  dplyr::mutate(name_in = str_c(file_accession, ".", format), 
                name_out = str_c("s_", biosample,
                                 ".e", experiment_group, 
                                 ".r", sample_group, 
                                 ".", read_group,
                                 "_", strand,
                                 ".", assembly, ".bw"), 
                name_out = str_replace_all(name_out, " ", "_")) %>% 
  dplyr::select(name_in, name_out) %>% 
  dplyr::mutate(rename_out = str_c("mv ", name_in, " ", name_out)) %$%
  rename_out %T>% 
  readr::write_lines(., path = "rename_bigwig.sh")
  