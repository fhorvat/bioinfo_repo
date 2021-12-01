### INFO: Rosalind - problem no. 004
### PROBLEM: Calculate rabbit population with mortal rabbits. 
###          Input: n - number of months, m - number of months rabbits live 
###          Each rabbit reproduces after 1 month and gives 1 offspring  
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
input_data <- read_lines(file.path(input_dir, "rosalind_005.txt"))

######################################################## MAIN CODE
# split input data to 2 vectors
input_data <- 
  stringr::str_split(input_data, pattern = " ") %>% 
  unlist() %>% 
  as.integer()

# rabbit(n, m) = rabbit(n - 1) + rabbit(n - 2) - rabbit(n - (m + 1))

fib <- function(n, m){
  
  # initialize vector
  x <- numeric(n)
  x[1:2] <- c(1, 1)
  
  # fill vector
  for(i in 3:n){
    if(i < (m + 1)){
      # normal Fibonacci sequence
      x[i] <- x[i - 1] + x[i - 2]
    }else{
      if(i == (m + 1)){ 
        # special case: i - (m + 1) < 1, add missing value
        x[i] <- x[i - 1] + x[i - 2] - 1
      }else{ 
        # i - (m + 1) >= 1
        x[i] <- x[i - 1] + x[i - 2] - x[i - (m + 1)]
      }
    }
  }
  
  return(x[n])
  
}

fib(input_data[1], input_data[2])
