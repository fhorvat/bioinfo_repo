### INFO: Rosalind - problem no. 004
### PROBLEM: calculate rabbit population. Input: n - number of months, k - number of rabbit offspring. Each rabbit reproduces after 1 month 
### DATE: 07. 06. 2017.  
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
input_data <- read_lines(file.path(input_dir, "rosalind_004.txt"))

######################################################## MAIN CODE
# split input data to 2 vectors
input_data <- 
  stringr::str_split(input_data, pattern = " ") %>% 
  unlist() %>% 
  as.integer()

# Fibbonaci rabbits
fib <- function(n, k) {
  
  # initialize vector
  x <- numeric(n)
  x[1:2] <- c(1, 1)
  
  # fill vector
  for(i in 3:n){
    x[i] = x[i - 1] + x[i - 2] * k 
  } 
  
  return(x[n])
}

fib(n = input_data[1], k = input_data[2])
