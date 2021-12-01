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
library(ggplot2)
library(scales)

################################################################################### SOURCE FILES

################################################################################### FUNCTIONS

################################################################################### PATH VARIABLES
# set in and out path
inpath <- getwd()

################################################################################### TABLES

################################################################################### MAIN CODE
### part 1 
readr::read_lines(file = "day01_input.txt") %>% 
  stringr::str_split(string = ., pattern = "") %>%
  unlist(.) %>% 
  tibble::tibble(input = .) %>%
  dplyr::mutate(input_lead = dplyr::lead(input), 
                input_lead = replace(input_lead, is.na(input_lead), input[1])) %>%
  dplyr::filter(input == input_lead) %>%
  dplyr::summarise(sum_input = sum(as.integer(input))) %$% 
  sum_input

### part 2
input_df <- 
  readr::read_lines(file = "day01_input.txt") %>% 
  stringr::str_split(string = ., pattern = "") %>%
  unlist(.) %>% 
  tibble::tibble(input = .) %>%
  dplyr::mutate(input_lead = dplyr::lead(input, n = (length(input) / 2)))

# replace lower half of lead column with upper half of original column
input_df[((nrow(input_df) / 2) + 1):nrow(input_df), 2] <- input_df[1:(nrow(input_df) / 2), 1]

# get answer
input_df %>%
  dplyr::filter(input == input_lead) %>%
  dplyr::summarise(sum_input = sum(as.integer(input))) %$% 
  sum_input