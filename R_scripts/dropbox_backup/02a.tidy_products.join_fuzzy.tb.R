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
## get common product manufactors in central and webshop tables
# intersect of names
product_manufact_shared <- intersect(webshop_tb$product_manufact, central_tb$product_manufact)


## set rules for string manipulation
# remove common words
common_words_remove <- "gratis|1\\+1|1 \\+ 1|1\\+ 1"

# replace untidy quantities
quantites_replace <- c("([0-9]{1,4}) (ml|g|mg|kom)" = "\\1\\2", "([0-9]{1,4})(ml|g|mg|kom)" = "\\1\\2 ",
                       "([a-z]+)([0-9]{1,4}(ml|g|mg))" = "\\1 \\2")

# replace some symbols (+ . - ... ) with whitespace
symbols_replace <- c("\\+|\\.|-|\\(|\\)|/" = " ")

# replace words with short universal version
short_words_replace <- c("djecj\\w* " = "djecj ", "vitamin\\w* " = "vitamin ", "zastitn\\w* " = "zastit ", 
                         "hydra\\w* " = "hidra ", "mineral*\\w " = "mineral ", "intenziv\\w* " = "intenziv ", 
                         "intim\\w* " = "intim ", "micelar\\w* " = "micelar ", "obnavljajuc\\w *" = "obnavljajuc ", 
                         "pjenus\\w* " = "pjenus ", "emolijent\\w* " = "emolijent ")

# replace vitamin with joined vitamin_x
vitamin_replace <- c("(vitamin) (a |b[0-9]+ |b |c |d |e |k |d3 |k2 |k1 )" = "\\1_\\2", 
                     "( a| b[0-9]+| b| c| d| e| k| d3| k2 | k1) (vitamin)" = "\\1_\\2")

# replace common words in webshop products
webshop_common_words_replace <- c("spray" = "sprej", "stick" = "stik")

# replace common words in central product
central_common_words_replace <- c(" kr |^kr " = " krema ", " tbl |^tbl | tbl$" = " tableta ", " caps |^caps " = " kapsule ", " kaps |^kaps " = " kapsule ", 
                                  " pr " = " protiv ", " mic ot |^mic ot " = " micelarna otopina ", " mic v |^mic v " = " micelarna voda ", 
                                  " deo |^deo " = " deodorant ", "^hyal | hyal " = " hyaluron ", "^los | los " = " losion ", " s/k " = " za suhu kozu ",
                                  " sum |^sum " = " sumece ", " cis |^cis | cisc |^cisc " = " ciscenje ", " mj mas " = " mjesovite masne ", 
                                  " spray |^spray " = " sprej ", "^samp | samp " = " sampon ", "^dj | dj " = " djecj ", " vit |^vit " = "vitamin", 
                                  " zast |^zast " = " zastit ", "hydra" = "hidra", " sun |^sun " = " suncanje ", " exp |^exp " = " express ",
                                  " int |^int " = " intenziv ", " obn |^obn " = " obnavljajuc ", "stick | ^stick" = " stik ", 
                                  " con emo " = " control emolijent ", 
                                  "gratis| grt$" = " ")


### tidy tables
# tidy webshop products
webshop_tidy <- 
  webshop_tb %>% 
  dplyr::select(product_id, product_manufact, product_title_tidy, webshop_price) %>% 
  dplyr::filter(!is.na(product_manufact)) %>% 
  dplyr::filter(product_manufact %in% product_manufact_shared) %>% 
  dplyr::mutate(product_title_tidy = 
                  product_title_tidy %>% 
                  str_remove_all(., common_words_remove) %>% 
                  str_replace_all(., c(quantites_replace, symbols_replace, webshop_common_words_replace, short_words_replace, vitamin_replace)) %>% 
                  str_squish(.)) %>% 
  dplyr::mutate(product_title_tidy =
                  product_title_tidy %>%
                  str_split(., " ") %>%
                  purrr::map(., sort) %>%
                  purrr::map(., str_c, collapse = " ") %>%
                  unlist(.))

# tidy central products
central_tidy <- 
  central_tb %>% 
  dplyr::select(product_title_tidy, product_manufact, central_price) %>% 
  dplyr::filter(product_manufact %in% product_manufact_shared) %>% 
  dplyr::mutate(product_title_tidy = 
                  product_title_tidy %>% 
                  str_replace_all(., c(quantites_replace, symbols_replace, central_common_words_replace, short_words_replace, vitamin_replace)) %>% 
                  str_squish(.)) %>% 
  dplyr::mutate(product_title_tidy =
                  product_title_tidy %>%
                  str_split(., " ") %>%
                  purrr::map(., sort) %>%
                  purrr::map(., str_c, collapse = " ") %>%
                  unlist(.))


### get final joined table and write it
# join using Jaccard distance
prices_tb <- 
  fuzzyjoin::stringdist_left_join(webshop_tidy[webshop_tidy$product_manufact == "la roche-posay", ], 
                                  central_tidy[central_tidy$product_manufact == "la roche-posay", ], 
                                  by = c("product_manufact", "product_title_tidy"), 
                                  method = "jaccard", 
                                  max_dist = 0.3) %>% 
  dplyr::select(product_id, product_manufact = product_manufact.x, 
                webshop_product = product_title_tidy.x, central_product = product_title_tidy.y,
                webshop_price, central_price) %>% 
  group_by(product_manufact, webshop_product) %>% 
  slice(1) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(same_prices = (round(webshop_price, 2) == round(central_price, 2)))

# count most common words
# central_tidy$product_title_tidy %>% str_c(., collapse = " ") %>% str_split(., " ") %>% table %>% as_tibble() %>% arrange(desc(n)) %>% asdf() %>% .[30:50, ]
# 
# central_tidy %>% dplyr::filter(str_detect(product_title_tidy, " rt |^rt "))
# 
# webshop_tidy %>% dplyr::filter(product_manufact == "apivita")
# 
# webshop_tidy %>% dplyr::filter(str_detect(product_title_tidy, "mjesovitu"))
# 
# webshop_tb %>% dplyr::filter(product_id %in% c(966, 974))



# # write table
# readr::write_csv(prices_tb, file.path(outpath, "prices_joined.tb.csv"))

