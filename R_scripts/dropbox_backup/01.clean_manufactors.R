### INFO: 
### DATE: Thu Jul 11 21:41:53 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/Tea")

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
outpath <- file.path(getwd(), "tidy")

# webshop prices path
webshop_path <- file.path(inpath, "Abeceda.xlsx")

# central prices path
central_path <- file.path(inpath, "central_prices.all.tsv")

######################################################## READ DATA
# read webshop prices
webshop <- openxlsx::read.xlsx(webshop_path)

# read central prices
central <- readr::read_delim(central_path, delim = "\t")

######################################################## MAIN CODE
### prepare some strings
# set manufactors with space in the name
webshop_manufact_with_space <- 
  c("dr theiss", "Dr. Brandt", "La Roche-Posay", "Natural Wealth", "Pharma Theiss", "Microlife", "Better You", 
    "Child Life") %>% 
  tolower(.)


#### clean webshop manufactors
# manufactors replacement strings
webshop_replace<- c("avene" = "avene", "chlidlife" = "childlife", "skinfnity" = "skinfinity",
                    "sforsin-" = "sforsin ", "^ricerfarma gengigel" = "gengigel ricerfarma", 
                    "è|æ" = "c", "š" = "s", "ž" = "z", "é" = "e", "â" = "a", "®" = "")

# create original table with product ID's
webshop_original <- 
  webshop %>% 
  as_tibble(.) %>% 
  dplyr::select(product_title = Title, webshop_price = Variant.Price) %>% 
  dplyr::mutate(product_id = 1:n())

# tidy webshop table
webshop_tidy <- 
  webshop_original %>%
  dplyr::mutate(product_title_tidy = 
                  product_title %>% 
                  tolower(.) %>%
                  str_replace_all(., webshop_replace)) %>% 
  dplyr::mutate(product_manufact = str_remove(product_title_tidy, " .*"),
                product_manufact_with_space = str_extract(product_title_tidy, str_c(webshop_manufact_with_space, collapse = "|")), 
                product_manufact = ifelse(is.na(product_manufact_with_space), product_manufact, product_manufact_with_space)) %>% 
  dplyr::mutate(product_title_tidy = str_remove(product_title_tidy, product_manufact)) %>% 
  dplyr::select(product_id, product_manufact, product_title_tidy) %>% 
  dplyr::filter(!str_detect(product_manufact, "\\.test"))

# join with original table and save
webshop_final <- 
  webshop_original %>% 
  dplyr::left_join(., webshop_tidy, by = "product_id") %T>%
  readr::write_csv(., path = file.path(outpath, "webshop_prices.tidy.csv"))


### clean central manufactors
# manufactors replacment strings
central_replace <- c("^d/" = "dietpharm ", "^av" = "avene", "^better" = "better you", 
                     "^biod." = "bioderma", "^isla-" = "isla ", "^nw" = "natural wealth", 
                     "^sf" = "skinfinity", "^hel" = "heliocare", "^lrp" = "la roche-posay", 
                     "^ptc" = "pharma theiss", "^np sforsin" = "sforsin np", "^child life" = "childlife", 
                     "è|æ" = "c", "š" = "s", "ž" = "z", "é" = "e", "â" = "a", "®" = "")

# create original table with product ID's
central_original <- 
  central %>% 
  as_tibble(.) %>% 
  dplyr::rename(central_price = price) %>% 
  dplyr::mutate(product_id = 1:n()) 

# tidy central table
central_tidy <- 
  central_original %>% 
  dplyr::mutate(product_title_tidy = 
                  product_title %>% 
                  tolower(.) %>%
                  str_replace_all(., central_replace)) %>% 
  dplyr::mutate(product_manufact = str_remove(product_title_tidy, " .*"), 
                product_manufact_with_space = str_extract(product_title_tidy, str_c(webshop_manufact_with_space, collapse = "|")), 
                product_manufact = ifelse(is.na(product_manufact_with_space), product_manufact, product_manufact_with_space), 
                product_manufact = ifelse(product_manufact %in% c("toplomjer", "cijev", "inhalator", "maska", "tlakomjer"), manufact_in, product_manufact)) %>% 
  dplyr::mutate(product_title_tidy = str_remove(product_title_tidy, product_manufact)) %>% 
  dplyr::select(product_id, product_manufact, product_title_tidy)

# join with original table and save
central_final <-
  central_original %>%
  dplyr::left_join(., central_tidy, by = "product_id") %>%
  dplyr::select(-manufact_in) %T>%
  readr::write_csv(., path = file.path(outpath, "central_prices.tidy.csv"))

