### INFO: 
### DATE: Thu Oct 10 17:04:49 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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
### IN AND OUT
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### COMMAND LINE ARGUMENTS
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)
dataset_path <- args$dataset_path
documentation_path <- args$documentation_path
feature_coordinates <- args$feature_coordinates


# coverage table path
coverage_path <- list.files(inpath, pattern = ".*\\.coverage\\.csv", full.names = T)

# RPM table path
rpm_path <- list.files(inpath, pattern = ".*\\.RPM\\.csv", full.names = T)

# original table path
line1_path <- list.files(documentation_path, pattern = "2019_RNAi_piRNA_file3_L1_table\\.csv|LINE1\\.5K_to_7k\\.ORFs.annotated_exons_and_L1Base\\.csv", full.names = T)

######################################################## READ DATA
# read coverage table
coverage_tb <- readr::read_csv(coverage_path)

# read RPM table
rpm_tb <- readr::read_csv(rpm_path)

# read LINE1 table
line1_tb <-
  readr::read_csv(line1_path) %>%
  dplyr::select_at(vars(-starts_with("X")))

######################################################## MAIN CODE
# join RPM and coverage
coverage_rpm <- 
  coverage_tb %>% 
  dplyr::left_join(., rpm_tb, by = "rmsk_id") %>% 
  dplyr::select_at(vars(rmsk_id, sort(colnames(.)[-1]))) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

# join with original table, save
if(basename(getwd()) == "L1Base"){
  
  # final table
  final_tb <- dplyr::left_join(line1_tb, coverage_rpm, by = c("L1base_id" = "rmsk_id"))
    
}else{
  
  if(basename(getwd()) == "rmsk"){
    
    # final table
    final_tb <- dplyr::left_join(line1_tb %>% dplyr::mutate(rmsk_id = as.character(rmsk_id)), coverage_rpm, by = "rmsk_id")
    
  }else{
    
    # warning 
    stop("You screwed something up mate! Check your paths and scripts!")
    
  }
  
}

# write table
readr::write_csv(final_tb, file.path(outpath, line1_path %>% basename(.) %>% str_replace(., "\\.csv", str_c(".", experiment, ".coverage_and_RPM.csv"))))
