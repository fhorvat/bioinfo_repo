### INFO: 
### DATE: Fri Jul 12 13:57:51 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/Tea/tidy")

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

library(fuzzyjoin)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# webshop prices path
webshop_path <- file.path(inpath, "webshop_prices.tidy.csv")

# central prices path
central_path <- file.path(inpath, "central_prices.tidy.csv")

######################################################## READ DATA
# read webshop prices
webshop_tb <- readr::read_csv(file = webshop_path)

# read central prices
central_tb <- readr::read_csv(central_path)

######################################################## MAIN CODE
# split webshop prices
webshop_list <- 
  webshop_tb %>% 
  dplyr::select(-product_title) %>% 
  dplyr::filter(!is.na(product_manufact)) %>% 
  split(., .$product_manufact)

# split central prices
central_list <- 
  central_tb %>% 
  dplyr::select(-product_title) %>% 
  split(., .$product_manufact)

# intersect of names
product_manufact_shared <- intersect(names(webshop_list), names(central_list))

# filter webshop and central list to include only shared manufacturer names
webshop_list <- webshop_list[names(webshop_list) %in% product_manufact_shared]
central_list <- central_list[names(central_list) %in% product_manufact_shared]


### join webshop and central tables
prices_tb <- 
  purrr::map(names(webshop_list), function(name){
    
    # get webshop prices for one manufact.
    webshop_filt <- 
      webshop_list[[name]] %>% 
      dplyr::select(product_id, product_manufact, product_title_tidy, webshop_price)
    
    # get central prices for one manufact.
    central_filt <- 
      central_list[[name]] %>% 
      dplyr::select(product_title_tidy, product_manufact, central_price)
    
    # join using Jaccard distance
    fuzzyjoin::stringdist_left_join(webshop_filt, central_filt, by = c("product_manufact", "product_title_tidy"), method = "jaccard", max_dist = 0.3) %>% 
      dplyr::mutate(product_manufact = name) %>% 
      dplyr::select(product_id, product_manufact, webshop_product = product_title_tidy.x, webshop_price, central_product = product_title_tidy.y, central_price) %>% 
      group_by(webshop_product) %>% 
      slice(1)

  }) %>% 
  dplyr::bind_rows(.)

# write
readr::write_csv(prices_tb, file.path(outpath, "prices_joined.list.csv"))



