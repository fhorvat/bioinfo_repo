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
input <- readr::read_lines(file = "day08_input.txt")

# rearrange input
exp_df <- 
  tibble(exp = input) %>% 
  tidyr::separate(col = exp, into = c("comm", "eval"), sep = "if") %>% 
  dplyr::select(2, 1) %>% 
  dplyr::mutate(comm = str_replace(comm, "inc", "+"),
                comm = str_replace(comm, "dec", "-")) %>% 
  dplyr::mutate(comm = str_c(str_extract(comm, "[[:alpha:]]+"), " <- ", comm), 
                var1 = str_extract(eval, "[[:alpha:]]+"), 
                var2 = str_extract(comm, "[[:alpha:]]+"))

# initialize all variables
variables_int <- unique(c(exp_df$var1, exp_df$var2))

for(x in 1:length(variables_int)){
  assign(variables_int[x], 0)
}

# go through data.frame and run commmands
for(x in 1:nrow(exp_df)){
  
  if(eval(parse(text = exp_df$eval[x]))){
    eval(parse(text = exp_df$comm[x]))
  }
  
}

# get largest variable 
sapply(variables_int, get, envir = .GlobalEnv) %>% 
  max


## part 2
# initialize all variables again
for(x in 1:length(variables_int)){
  assign(variables_int[x], 0)
}

# go through data.frame and run commmands, check in every step which is largest variable
max_value <- NULL

for(x in 1:nrow(exp_df)){
  
  if(eval(parse(text = exp_df$eval[x]))){
    eval(parse(text = exp_df$comm[x]))
  }
  
  max_value <- c(max_value, sapply(variables_int, get, envir = .GlobalEnv) %>% max)
  
}

max(max_value)
