### INFO: Rosalind - problem no. 002
### PROBLEM: transform DNA string to RNA string (replace T with U)
### DATE: 06. 06. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/rosalind")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)

######################################################## PATH VARIABLES
input_dir <- "C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/rosalind/input"

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
input_data <- read_lines(file.path(input_dir, "rosalind_002.txt"))

######################################################## MAIN CODE
# split string to vectors of length 1, count unique characters
input_data %>% 
  stringr::str_replace_all(., pattern = "T", replacement = "U")

