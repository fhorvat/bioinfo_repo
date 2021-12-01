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
  readr::read_delim(file = "day06_input.txt", delim = "\t", col_names = F) %>% 
  as.integer(.)

state_all <- "x"

while(!any(duplicated(state_all))){
  
  max_pos <- which.max(input)
  blocks <- input[max_pos]
  input[max_pos] <- 0
  
  n_all <- blocks %/% length(input)
  input <- input + n_all
  blocks <- blocks %% length(input)
  
  if(blocks > 0){
    
    if(max_pos < length(input)){
      input[sort(c((max_pos + 1):length(input), 1:max_pos)[1:blocks])] <-  input[sort(c((max_pos + 1):length(input), 1:max_pos)[1:blocks])] + 1 
    }else{
      input[c(1:max_pos)[1:blocks]] <-  input[c(1:max_pos)[1:blocks]] + 1 
    }
    
  }
  
  state <- str_c(input, collapse = " ")
  state_all <- c(state_all, state)
  
}

length(state_all) - 1

### part 2
cycle <- which(state_all == state_all[duplicated(state_all)])
cycle[2] - cycle[1]
