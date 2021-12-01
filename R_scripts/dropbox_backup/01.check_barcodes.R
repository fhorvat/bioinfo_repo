### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Raw/Links")

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

# barcodes count path
barcodes_path <- list.files(inpath, "\\.barcodes\\.txt", full.names = T)

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Documentation"

# sample table path
sample_table_path <- list.files(documentation_path, ".*\\.sampleTable\\.csv", full.names = T)

######################################################## READ DATA
# read most common barcodes for each sample
barcodes_tb <- purrr::map(barcodes_path, function(path){
  
  readr::read_delim(path, col_names = c("count", "barcode_seq"), delim = " ") %>% 
    dplyr::mutate(count = count %>% str_trim(.) %>% as.numeric(.)) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.barcodes\\.txt"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(sample_id, barcode_seq, count)

# read sample table
sampleTable <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
# get the most 2 common barcodes per sample
barcodes_top_2 <- 
  barcodes_tb %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::top_n(., n = 2, wt = count) %>% 
  dplyr::ungroup(.)

# get the most common barcode
barcodes_top <- 
  barcodes_tb %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::filter(count == max(count)) %>% 
  dplyr::ungroup(.)

# join sample table with most common barcode, check if they match
sample_tb_barcode <- 
  sampleTable %>% 
  dplyr::left_join(., barcodes_top, by = "sample_id") %>%
  dplyr::filter(barcode != barcode_seq)

# get percentage of other barcodes in each sample
barcodes_other <-
  barcodes_tb %>% 
  # dplyr::anti_join(., barcodes_top, by = c("sample_id", "barcode_seq")) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(count_other = sum(count[count != max(count)])) %>% 
  dplyr::filter(count == max(count)) %>% 
  dplyr::mutate(percentage_other = 100 * (count_other / (count + count_other)))
