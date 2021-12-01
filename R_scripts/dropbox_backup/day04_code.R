### INFO: Advent of Code 2018, day 01
### DATE: Sat Dec 01 13:44:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/test/adventOfCode/2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(lubridate)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# read data
input <- readr::read_lines(file = file.path(inpath, "day04_input.txt"))

######################################################## MAIN CODE
# split input to table
input_tb <- 
  tibble(raw_input = input) %>% 
  dplyr::mutate(date_time = str_extract(raw_input, "\\[.*\\]") %>% str_remove_all(., "\\[|\\]"), 
                description = str_remove(raw_input, "\\[.*\\] "), 
                date = str_remove(date_time, " .*")) %>% 
  dplyr::select(date_time, date, description) %>% 
  dplyr::mutate(date_time = ymd_hm(date_time),
                date = ymd(date)) %>% 
  dplyr::arrange(date_time)



