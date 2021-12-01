### INFO: Advent of code 2017, day01 (http://adventofcode.com/2017/day/1)
### DATE: 18. 12. 2017 
### AUTHOR: Filip Horvat

################################################################################### SCRIPT PARAMS
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/test/adventOfCode/2017")

################################################################################### LIBRARIES
# data shaping
library(magrittr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(purrrlyr)

################################################################################### SOURCE FILES

################################################################################### FUNCTIONS

################################################################################### PATH VARIABLES
# set in and out path
inpath <- getwd()

################################################################################### TABLES

################################################################################### MAIN CODE
### part 1 
readr::read_delim(file = "day02_input.txt", delim = "\t", col_names = F) %>% 
  purrrlyr::by_row(., max, .collate = "cols", .to = "max") %>% 
  purrrlyr::by_row(., min, .collate = "cols", .to = "min") %>% 
  dplyr::mutate(diff = max - min) %>% 
  dplyr::summarise(checksum = sum(diff)) %$%
  checksum

# dt <- 
#   readr::read_delim(file = "day02_input.txt", delim = "\t", col_names = F) %>% 
#   data.table::as.data.table(.)
# dt[, `:=` (MIN = rowMins(as.matrix(.SD), na.rm = T),
#            MAX = rowMaxs(as.matrix(.SD), na.rm = T))]

### part 2
# find pair of numbers in each row which are evenly divisible 
rows <- readr::read_delim(file = "day02_input.txt", delim = "\t", col_names = F)

lapply(X = 1:nrow(rows), FUN = function(X){
  as.integer(rows[X, ]) %>% 
    expand.grid(., .) %>% 
    dplyr::filter(Var1 != Var2) %>% 
    dplyr::mutate(x_y_div = Var1 / Var2, 
                  x_y_modulo = Var1 %% Var2) %>% 
    dplyr::filter(x_y_modulo == 0) %$%
    x_y_div
}) %>% 
  unlist(.) %>% 
  sum(.)



