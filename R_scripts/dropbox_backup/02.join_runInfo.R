### INFO: 
### DATE: Wed May 02 13:07:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/sequencing_info/RNAseq")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA
# read SRA tables - overview of experiments and details about samples 
sra_overview <- readr::read_csv(file = file.path(inpath, "SRA2018_overview.csv"))
sra_details <- readr::read_csv(file = file.path(inpath, "SRA2018_details.csv"))

######################################################## MAIN CODE
# join details with overview
sra_all <- 
  left_join(sra_details, sra_overview %>% dplyr::select(which(!colnames(.) %in% colnames(sra_details))), 
            by = c("experiment_id" = "reference")) %>% 
  dplyr::mutate(library_strand = ifelse(is.na(library_strand), library_strand_specific, library_strand)) %>% 
  dplyr::select(experiment_id:library_scope, library_selection, library_layout, library_strand, read_length, Instrument, everything()) %>% 
  dplyr::select(-c(assay, organism, replicate_number, library_platform, raw_path, mapped_path, 
                   mapped_genome, X25, library_strand_specific))

# save
readr::write_csv(x = sra_all, path = file.path(outpath, "SRA2018.RNAseq.csv"))

# get all cells/timepoints per experiment
sra_all %>% 
  dplyr::group_by(experiment_id) %>% 
  dplyr::distinct(cell, time_point, .keep_all = T) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(time_point = replace(time_point, is.na(time_point), "")) %>% 
  dplyr::mutate(cell_timepoint = str_c(cell, " ", time_point) %>% stringr::str_trim(.)) %>% 
  dplyr::group_by(experiment_id) %>% 
  dplyr::summarise(cell_timepoint = str_c(unique(cell_timepoint), collapse = '|'), 
                   genotype = str_c(unique(genotype), collapse = '|'), 
                   treatment = str_c(unique(treatment), collapse = "|"), 
                   library_scope = str_c(unique(library_scope), collapse = "|")) %T>% 
  readr::write_csv(., path = file.path(outpath, "SRA2018.RNAseq.collapsed.csv"))

  
