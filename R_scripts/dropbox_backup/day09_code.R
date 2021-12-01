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
  readr::read_lines(file = "day09_input.txt") %>% 
  str_split(., pattern = "", simplify = F) %>% 
  unlist(.)

input <- "{{<ab>},{<ab>},{<ab>},{<ab>}}" %>% 
  str_split(., pattern = "", simplify = F) %>% 
  unlist(.)

group <- 0

for(x in 1:length(input)){
  
  if(!ignore.next){
    
    if(input[x] == "!"){
      ignore.next <- T
    }
    
    if(input[x] == "{"){
      group <- group + 1
    }
    
    if(input[x] == "}"){
      group <- group - 1
    }
    
  }
  
  
}

