### INFO: 
### DATE: Wed May 02 13:07:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/sequencing_info/RNAseq/runInfo")

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

# runInfo paths
runinfo_path <- list.files(path = inpath, pattern = "*.runInfo.txt", full.names = T)

######################################################## READ DATA
# read all runInfo tables
runinfo_df <- 
  lapply(runinfo_path, function(x){
    
    # read and filter runInfo, add info about experiment
    readr::read_delim(file = x, delim = "\t", col_types = cols(.default = "c")) %>% 
      dplyr::mutate(experiment_id = basename(x) %>% stringr::str_remove(., ".runInfo.txt"))
    
  }) %>% 
  dplyr::bind_rows(.)
  
######################################################## MAIN CODE
# save whole table
readr::write_csv(runinfo_df, path = file.path(outpath, "all_runInfo.csv"))

runinfo_df %>% 
  dplyr::filter(!is.na(genetic_background)) %>% 
  dplyr::select(which(colSums(apply(., 2, is.na)) == 0)) %>%
  asdf

runinfo_df %>% 
  dplyr::filter(experiment_id == "Ramskold_2012_NatBiotechnol_GSE38495") %>% 
  dplyr::select(which(colSums(apply(., 2, is.na)) == 0)) %>%
  asdf
