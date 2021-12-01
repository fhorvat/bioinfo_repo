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
input <- 
  readr::read_lines(file = "day05_input.txt") %>% 
  as.integer(.)
pos <- 1
step <- 0

while(pos <= length(input)){
  
  jump <- input[pos]
  input[pos] <- input[pos] + 1
  pos <- pos + jump
  step <- step + 1
  
}

print(step)

### part 2 
input <- 
  readr::read_lines(file = "day05_input.txt") %>% 
  as.integer(.)
pos <- 1
step <- 0

while(pos <= length(input)){
  
  jump <- input[pos]
  
  if(jump >= 3){
    input[pos] <- input[pos] - 1
  }else{
    input[pos] <- input[pos] + 1
  }

  pos <- pos + jump
  step <- step + 1
  
}

print(step)
