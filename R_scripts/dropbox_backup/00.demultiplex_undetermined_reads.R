### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq/Data/Documentation")

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

# raw files path
raw_path <- "/common/RAW/Svoboda/golden_hamster.testis.smallRNAseq.8.5dpp.2020/02_undetermined_reads"

# barcode lists path
barcode_path_list <- list.files(raw_path, "\\.barcodes\\.txt", full.names = T)

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw\\.csv", full.names = T)

######################################################## READ DATA
# read the barcode list
barcode_tb <- purrr::map(barcode_path_list, function(path){
  
  readr::read_delim(path, delim = " ", col_names = c("count", "barcode")) %>% 
    dplyr::mutate(count = count %>% as.numeric(.), 
                  lane = path %>% basename(.) %>% str_extract(., "L00[1-4]"))
  
}) %>% 
  bind_rows(.)

# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# get top 4 barcodes per lane
barcode_top4 <- 
  barcode_tb %>% 
  dplyr::group_by(lane) %>% 
  dplyr::slice_max(order_by = count, n = 4) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(barcode)

# check if those are all present in sample table
all(barcode_top4$barcode %in% sample_table_raw$barcode)

# save barcodes
sample_table_raw$barcode %>% 
  unique(.) %>% 
  readr::write_lines(., path = file.path(outpath, "barcodes_clean.txt"))

