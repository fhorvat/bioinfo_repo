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
input <- 289326

## create spiral matrix
# initialize
mat <- matrix(c(1, 2), nrow = 1)
start <- 3
seq_length <- 2

# expand matrix
while(!any(mat == input)){
  
  mat <- rbind(rev(seq.int(from = start, by = 1, length.out = seq_length)), mat)
  start <- max(mat) + 1
  
  mat <- cbind(seq.int(from = start, by = 1, length.out = seq_length), mat)
  start <- max(mat) + 1
  seq_length <- seq_length + 1 
  
  mat <- rbind(mat, seq.int(from = start, by = 1, length.out = seq_length))
  start <- max(mat) + 1
  
  mat <- cbind(mat, rev(seq.int(from = start, by = 1, length.out = seq_length)))
  start <- max(mat) + 1
  seq_length <- seq_length + 1 
  
}

# get distance from input to start position
pos_mat <- rbind(which(mat == input, arr.ind = T), 
                 which(mat == 1, arr.ind = T))
row_dist <- abs(pos_mat[1, 1] - pos_mat[2, 1])
col_dist <- abs(pos_mat[1, 2] - pos_mat[2, 2])
row_dist + col_dist


### part 2
mat <- matrix(c(1, 1), nrow = 1) 
n_row <- 1
n_col <- 1

## 1.
repeat{
  
  mat <- rbind(0, mat)
  
  for(x in ncol(mat):1){
    
    if(mat[n_row, n_col] > input){
      break   
    }
    
    n_row <- 1
    all_rows <- (n_row - 1):(n_row + 1)
    all_rows <- all_rows[(all_rows <= nrow(mat)) & (all_rows > 0)]
    
    n_col <- x
    all_cols <- (n_col - 1):(n_col + 1)
    all_cols <- all_cols[(all_cols <= ncol(mat)) & (all_cols > 0)]
    
    mat[n_row, n_col] <- sum(mat[all_rows, all_cols])
    
  }
  
  if(mat[n_row, n_col] > input){
    output <- mat[n_row, n_col]
    break   
  }
  
  ## 2.
  mat <- cbind(0, mat)
  
  for(x in 1:nrow(mat)){
    
    if(mat[n_row, n_col] > input){
      break   
    }
    
    n_row <- x
    all_rows <- (n_row - 1):(n_row + 1)
    all_rows <- all_rows[(all_rows <= nrow(mat)) & (all_rows > 0)]
    
    n_col <- 1
    all_cols <- (n_col - 1):(n_col + 1)
    all_cols <- all_cols[(all_cols <= ncol(mat)) & (all_cols > 0)]
    
    mat[n_row, n_col] <- sum(mat[all_rows, all_cols])
    
  }
  
  if(mat[n_row, n_col] > input){
    output <- mat[n_row, n_col]
    break   
  }
  
  ## 3.
  mat <- rbind(mat, 0)
  
  for(x in 1:ncol(mat)){
    
    if(mat[n_row, n_col] > input){
      break   
    }
    
    n_row <- nrow(mat)
    all_rows <- (n_row - 1):(n_row + 1)
    all_rows <- all_rows[(all_rows <= nrow(mat)) & (all_rows > 0)]
    
    n_col <- x
    all_cols <- (n_col - 1):(n_col + 1)
    all_cols <- all_cols[(all_cols <= ncol(mat)) & (all_cols > 0)]
    
    mat[n_row, n_col] <- sum(mat[all_rows, all_cols])
  
  }
  
  if(mat[n_row, n_col] > input){
    output <- mat[n_row, n_col]
    break   
  }
  
  ## 4.
  mat <- cbind(mat, 0)
  
  for(x in nrow(mat):1){
    
    if(mat[n_row, n_col] > input){
      break   
    }
    
    n_row <- x
    all_rows <- (n_row - 1):(n_row + 1)
    all_rows <- all_rows[(all_rows <= nrow(mat)) & (all_rows > 0)]
    
    n_col <- ncol(mat)
    all_cols <- (n_col - 1):(n_col + 1)
    all_cols <- all_cols[(all_cols <= ncol(mat)) & (all_cols > 0)]
    
    mat[n_row, n_col] <- sum(mat[all_rows, all_cols])
    
  }
  
  if(mat[n_row, n_col] > input){
    output <- mat[n_row, n_col]
    break   
  }
  
}

print(output)
