### INFO: 
### DATE: Thu Jul 11 21:41:53 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/Tea/central_prices")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# central prices path
central_path <- list.files(path = inpath, pattern = ".*\\.csv", full.names = T)

######################################################## READ DATA
# read all .csv files, tidy, join together
central_tb <- purrr::map(central_path, function(path){
  
  # read and tidy one file
  suppressMessages(central_tidy <- 
                     readr::read_delim(path, delim = "\t", skip = 12, col_names = "tmp", skip_empty_rows = T) %>% 
                     dplyr::mutate(tmp = str_trim(tmp)) %>% 
                     tidyr::separate(tmp, into = c("product_title", "price"), sep = " {7,}") %>% 
                     dplyr::mutate(manufact_in = basename(path) %>% str_remove(., "\\.csv"))) 
  
}) %>% 
  bind_rows(.) %>% 
  dplyr::filter(!(is.na(price)), 
                !(str_detect(price, "Datum|Stranica"))) %>% 
  dplyr::mutate(price = str_replace(price, ",", "."))

######################################################## MAIN CODE
# save table
readr::write_delim(central_tb, path = file.path(inpath, "..", "central_prices.all.tsv"), delim = "\t")



