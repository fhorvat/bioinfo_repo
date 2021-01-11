### INFO: 
### DATE: Wed Sep 30 13:28:24 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
fastq_name <- args$fastq_name

######################################################## READ DATA
# read first read in a pair
read_1_list <- readr::read_lines(file = file.path(inpath, str_c(fastq_name, ".PE_1.converted_reads.txt")))

# read second read in a pair
read_2_list <- readr::read_lines(file = file.path(inpath, str_c(fastq_name, ".PE_2.converted_reads.txt")))

######################################################## MAIN CODE
# intersect two lists and save to disk
intersect(read_1_list, read_2_list) %>% 
  readr::write_lines(., file.path(outpath, str_c(fastq_name, ".converted_reads.intersect.txt")))
