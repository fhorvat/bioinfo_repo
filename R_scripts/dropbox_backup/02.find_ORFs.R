### INFO: 
### DATE: Thu Nov 07 08:43:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/hackaton")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### read and parse fasta
read_fasta <- function(path){
  
  # read lines
  fasta_raw <- readr::read_lines(path)
  
  # add separator at the end of the header
  which_header <- which(str_detect(fasta_raw, ">"))
  fasta_raw[which_header] <- str_c(fasta_raw[which_header], "<")
  
  # concatenate sequences together
  fasta_cat <- str_c(fasta_raw, collapse = "")
  
  # split by ">"
  fasta_split <- 
    str_split(fasta_cat, pattern = ">") %>% 
    unlist(.) %>% 
    .[-1]
  
  # get and remove names
  fasta_names <- str_extract(fasta_split, ".*(?=\\<)")
  fasta_split <- str_remove(fasta_split, ".*\\<")
  
  # set names and return
  names(fasta_split) <- fasta_names
  
  # return
  return(fasta_split)
  
}

### reverse complement
revComp <- function(string, reverse = T){
  
  # split string
  string_split <- str_split(string, "") %>% unlist(.)
  
  # reverse
  if(reverse){
    string_out <- rev(string_split)
  }else{
    string_out <- string_split
  }
  
  # complement
  if(complement){
    string_out <- 
      factor(string_out, levels = c("A", "C", "G", "T"), labels = c("T", "G", "C", "A")) %>% 
      as.character(.)
  }
  
  # paste together
  string_out <- str_c(string_out, collapse = "")
  
  # return
  return(string_out)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# fasta path
fasta_path <- file.path(inpath, "L1_full_length.chr11.fasta")
  
######################################################## READ DATA
# read fasta
fasta_raw <- read_fasta(fasta_path)

######################################################## MAIN CODE

fasta_rc <- purrr::map(fasta_raw, revComp) %>% unlist(.)





