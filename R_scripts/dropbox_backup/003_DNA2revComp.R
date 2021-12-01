### INFO: Rosalind - problem no. 003
### PROBLEM: make a reverse complement of DNA string 
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
input_data <- read_lines(file.path(input_dir, "rosalind_003.txt"))

######################################################## MAIN CODE
# reverse complement
input_data %>% 
  stringr::str_split(., pattern = "") %>% # split string
  unlist(.) %>% 
  sapply(., switch, "A" = "T", "T" = "A", "G" = "C", "C" = "G") %>% # complement
  unname() %>% 
  stringr::str_c(., collapse = "") %>% # collapse string
  stringi::stri_reverse(.) # reverse
