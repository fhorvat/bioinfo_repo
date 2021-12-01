### INFO: reads .out files from repeatMasker
### DATE: 13. 02. 2018.
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/repeatMasker_hmmer")

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
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# list all .out tables 
out_tables <- list.files(path = file.path(outpath, "out"), pattern = "*.fasta.out", full.names = T)

######################################################## MAIN CODE
# clean rmsk.out table, write
for(out_table in out_tables){
  
  readr::read_table2(file = out_table, skip = 3, col_names = F) %>%
    dplyr::select(query = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, id = X15) %>%
    tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
    dplyr::mutate(strand = replace(strand, strand == "C", "*")) %>% 
    dplyr::filter(repClass == "LTR") %T>%
    readr::write_delim(., path = str_replace(out_table, ".fasta.out", ".clean.fasta.out") %>% basename(.)) %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise(query = unique(query), 
                     start = min(start),
                     end = max(end), 
                     repName = str_c(repName, collapse = "|"), 
                     repClass = str_c(unique(repClass), collapse = "|"), 
                     repFamily = str_c(unique(repFamily), collapse = "|")) %>% 
    dplyr::select(-id) %T>% 
    readr::write_delim(., path = str_replace(out_table, ".fasta.out", ".sum.clean.fasta.out") %>% basename(.))
  
}


