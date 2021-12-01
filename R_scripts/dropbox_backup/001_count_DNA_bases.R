### INFO: Rosalind - problem no. 001
### PROBLEM: count A, C, G, T bases in given DNA string 
### DATE: 05. 06. 2017.  
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
input_data <- read_lines(file.path(input_dir, "rosalind_001.txt"))

######################################################## MAIN CODE
# split string to vectors of length 1, count unique characters
input_data %>% 
  stringr::str_split(., pattern = "") %>% 
  table()

