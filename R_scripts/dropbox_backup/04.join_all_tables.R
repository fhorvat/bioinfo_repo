### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/library_composition")

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

# find all tables with RPM values
rpm_tb_path <- list.files(inpath, ".*\\.RPM\\.complete\\.csv", full.names = T)

readr::read_csv("LINE1s_subsets.testis.RNAseq.RPM.complete.csv")

######################################################## READ DATA
# read and join table
rpm_tb <- purrr::map(rpm_tb_path, function(path){
  
  # read table
  rpm <- readr::read_csv(path)
  
  # change subset for small RNA-seq (wrong name)
  rpm %<>% 
    dplyr::mutate(subset = replace(subset, subset == "LINE1_intact", "LINE1_FLIs"))
  
  # return
  return(rpm)
  
}) %>% 
  purrr::reduce(., left_join, by = "subset")

# save
readr::write_csv(rpm_tb, file.path(outpath, "LINE1s_subsets.all.RPM.csv"))

######################################################## MAIN CODE
