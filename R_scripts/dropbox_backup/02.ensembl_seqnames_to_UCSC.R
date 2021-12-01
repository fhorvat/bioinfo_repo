### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/test/mouse/mm10.GRCm38.p5.GCA_000001635.7")

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

# assembly report path
assembly_report_path <- list.files(path = inpath, pattern = ".*assembly_report.txt")

######################################################## READ DATA
# read assembly report
assembly_report <- readr::read_delim(file = assembly_report_path, delim = "\t", comment = "#", col_names = F)
  
######################################################## MAIN CODE
# add ensembl names to assembly report, write table
assembly_report %<>%
  magrittr::set_colnames(., c("sequence_name", "sequence_role", "assigned_molecule", "assigned_molecule_location_type", "genBank_accn", "relationship", 
                              "refSeq_accn", "assembly_unit", "sequence_length", "UCSC_name")) %>% 
  dplyr::mutate(ensembl_name = ifelse(sequence_role == "assembled-molecule", assigned_molecule, genBank_accn), 
                UCSC_name = replace(UCSC_name, assigned_molecule_location_type == "Mitochondrion", "chrM")) %>% 
  dplyr::select(ensembl_name, UCSC_name, sequence_role) %T>% 
  readr::write_delim(x = ., 
                     path = file.path(outpath, basename(assembly_report_path) %>% stringr::str_replace(., "_assembly_report.txt", ".ensembl2UCSC.txt")), 
                     delim = "\t")
  

