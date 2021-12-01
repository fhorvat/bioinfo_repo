### INFO: Advent of code 2017, day01 (http://adventofcode.com/2017/day/1)
### DATE: 18. 12. 2017 
### AUTHOR: Filip Horvat

################################################################################### SCRIPT PARAMS
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("C:/Users/Filip/Dropbox/Bioinfo/other_projects/test/adventOfCode/2017")

################################################################################### LIBRARIES
# data shaping
library(magrittr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)

################################################################################### SOURCE FILES

################################################################################### FUNCTIONS

################################################################################### PATH VARIABLES
# set in and out path
inpath <- getwd()

################################################################################### TABLES

################################################################################### MAIN CODE
### part 1 
readr::read_lines(file = "day04_input.txt") %>% 
  stringr::str_split(string = ., pattern = " ") %>% 
  lapply(X = ., function(X) !any(duplicated(X))) %>% 
  unlist(.) %>% 
  sum(.)

### part 2 
readr::read_lines(file = "day04_input.txt") %>% 
  stringr::str_split(string = ., pattern = " ") %>% 
  sapply(X = ., function(X){
    stringr::str_split(string = X, pattern = "") %>% 
      lapply(X = ., function(X){
        sort(X) %>% 
          str_c(., collapse = "")
      }) %>% 
      duplicated(.) %>% 
      any(.)
  }) %>% 
  not(.) %>% 
  sum(.)
