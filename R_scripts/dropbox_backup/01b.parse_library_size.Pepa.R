### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Documentation")

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

# RDA path 
rda_path <- file.path(inpath, "allPSg_objs.rda")

# sample table path
sample_tb_path <- file.path(inpath, "hamster_testis_Mov10l.20191008.sampleTable.csv")

######################################################## READ DATA
# load object 
load(rda_path, verbose = T)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# variable names
table_names <- ls()[str_detect(ls(), ".*\\.psg\\.dt")]

# get variable from environment
tb_sum <- purrr::map(table_names, function(tb){

  # get table, filter, sum rows, store in a tibble
  tb_sum <- 
    get(tb) %>% 
    tibble::as_tibble(.) %>% 
    # dplyr::filter(rLen >= 19, rLen <= 32) %>% 
    tidyr::pivot_longer(-rLen, names_to = "class", values_to = "count") %>% 
    # dplyr::group_by(rLen) %>% 
    dplyr::summarise(count_19to32nt.Pepa = sum(count)) %>% 
    dplyr::mutate(sample_name = tb %>% str_remove_all(., "HET_|KO_|WT_|\\.psg\\.dt") %>% str_replace(., "_", "-"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(sample_name, count_19to32nt.Pepa)
  
# add to sample table
sample_tb %<>%
  dplyr::left_join(., tb_sum, by = "sample_name")

# save 
readr::write_csv(sample_tb, sample_tb_path)

