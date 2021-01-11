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
fastq_path <- args$fastq_path

# get fastq name
fastq_name <- fastq_path %>% basename(.) %>% str_remove(., "\\.fastq")

######################################################## READ DATA

######################################################## MAIN CODE
# create a function which returns names of reads which pass conversion filter
getConvertedReadNames <- function(fastq_chunk, dummyVariable){
  
  # get DNAStringSet of reads seqeuce
  fastq_seq <- 
    fastq_chunk[seq(2, length(fastq_chunk), 4)] %>% 
    DNAStringSet(.)
  
  # nucleotide frequency in the table
  freq_tb <- 
    fastq_seq %>% 
    alphabetFrequency(., as.prob = T) %>% 
    .[, 1:4]
  
  # get reads which have < 5% C or G nucleotides
  converted_reads <- 
    fastq_chunk[seq(1, length(fastq_chunk), 4)][which(freq_tb[, 2] < 0.05 | freq_tb[, 3] < 0.05)] %>% 
    str_remove(., "^@")
  
  # save to file on disk
  readr::write_lines(x = converted_reads, path = file.path(outpath, str_c(fastq_name, "converted_reads", "txt", sep = ".")), append = T)
  
  # return
  return(NULL)

}

# run function on chunks of fastq file
readr::read_lines_chunked(fastq_path, SideEffectChunkCallback$new(getConvertedReadNames), chunk_size = 1000000)
