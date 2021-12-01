### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/library_composition")

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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# RPM histograms path
rpm_paths <- list.files(inpath, ".*read_length.histogram.csv", full.names = T)

######################################################## READ DATA
# read and summarize different read lengths
rpm_tb <- purrr::map(rpm_paths, function(path){
  
  # subset name
  subset_name <- 
    path %>% 
    basename(.) %>% 
    str_remove(., "\\.testis\\.small_RNAseq\\.RPM\\.read_length\\.histogram.csv") %>% 
    str_remove(., "^0[1-9]{1}\\.")
  
  # read table
  rpm <- 
    readr::read_csv(path) %>% 
    tidyr::pivot_longer(cols = -genotype_age, names_to = "width", values_to = "RPM") %>% 
    dplyr::mutate(width = width %>% str_remove(., "r\\.") %>% as.numeric(.)) %>% 
    dplyr::filter(width >= 24, width <= 31) %>% 
    dplyr::group_by(genotype_age) %>% 
    dplyr::summarise(RPM.24to31nt = sum(RPM)) %>%
    dplyr::mutate(subset = subset_name)
    
}) %>% 
  bind_rows(.) %>% 
  dplyr:::mutate(RPM.24to31nt = replace(RPM.24to31nt, is.na(RPM.24to31nt), 0)) %>% 
  tidyr::pivot_wider(., id_cols = subset, names_from = genotype_age, values_from = RPM.24to31nt, names_prefix = "RPM.24to31nt.") %>% 
  dplyr::mutate(rpm_tb = replace(rpm_tb, rpm_tb == "intact_guys", "IAP_FLIs")) %T>%
  readr::write_csv(., file.path(outpath, "IAPs_subsets.testis.small_RNAseq.RPM.complete.csv"))
  
######################################################## MAIN CODE
