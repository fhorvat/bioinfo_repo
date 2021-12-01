# --- Day 2: Inventory Management System ---
#   You stop falling through time, catch your breath, and check the screen on the device. "Destination reached. Current Year: 1518. Current Location: North Pole Utility Closet 83N10." You made it! Now, to find those anomalies.
# 
# Outside the utility closet, you hear footsteps and a voice. "...I'm not sure either. But now that so many people have chimneys, maybe he could sneak in that way?" Another voice responds, "Actually, we've been working on a new kind of suit that would let him fit through tight spaces like that. But, I heard that a few days ago, they lost the prototype fabric, the design plans, everything! Nobody on the team can even seem to remember important details of the project!"
# 
# "Wouldn't they have had enough fabric to fill several boxes in the warehouse? They'd be stored together, so the box IDs should be similar. Too bad it would take forever to search the warehouse for two similar box IDs..." They walk too far away to hear any more.
# 
# Late at night, you sneak to the warehouse - who knows what kinds of paradoxes you could cause if you were discovered - and use your fancy wrist device to quickly scan every box and produce a list of the likely candidates (your puzzle input).
# 
# To make sure you didn't miss any, you scan the likely candidate boxes again, counting the number that have an ID containing exactly two of any letter and then separately counting those with exactly three of any letter. You can multiply those two counts together to get a rudimentary checksum and compare it to what your device predicts.
# 
# For example, if you see the following box IDs:
# 
# abcdef contains no letters that appear exactly two or three times.
# bababc contains two a and three b, so it counts for both.
# abbcde contains two b, but no letter appears exactly three times.
# abcccd contains three c, but no letter appears exactly two times.
# aabcdd contains two a and two d, but it only counts once.
# abcdee contains two e.
# ababab contains three a and three b, but it only counts once.
# Of these box IDs, four of them contain a letter which appears exactly twice, and three of them contain a letter which appears exactly three times. Multiplying these together produces a checksum of 4 * 3 = 12.
# 
# What is the checksum for your list of box IDs?

### INFO: Advent of Code 2018, day 02
### DATE: Sat Dec 01 13:44:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/test/adventOfCode/2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# read data
input <- readr::read_lines(file = file.path(inpath, "day02_input.txt"))

######################################################## MAIN CODE
#### PART 1 ####
# count letters in each IDs
x <- 
  input %>% 
  str_split(., "") %>%
  purrr::map(., function(x){
    
    table(x) %>% 
      as.tibble %>% 
      dplyr::filter(n == 2 | n == 3) %$%
      n %>% 
      unique(.)
    
  }) %>% 
  unlist(.) %>% 
  table(.) %>% 
  as.tibble(.) %$%
  n 

# multiply to get checksum
`*`(x[1], x[2])


#### PART 2 ####
# get all combinations of strings
input_comb <- expand.grid(input, input, stringsAsFactors = F) 

# split
input1 <- str_split(input_comb$Var1, "")
input2 <- str_split(input_comb$Var2, "")

# in loop walk through string
counter <- 0
while(T){
  
  # increase counter
  counter <- counter + 1
 
  # compare inputs
  if(length(input1[[counter]][input1[[counter]] == input_2[[counter]]]) == 25){
    stop(str_c(input1[[counter]][input1[[counter]] == input_2[[counter]]], collapse = ""))
  }

}



