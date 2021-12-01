### INFO: 
### DATE: Wed Sep 30 13:28:24 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Raw/Links/sample_fastq")

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

# fastq path
fastq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Raw/Links/sample_fastq/sample.1.fastq"
  
# get fastq name
fastq_name <- fastq_path %>% basename(.) %>% str_remove(., "\\.fastq")

######################################################## READ DATA
# read fastq
fastq_file <- readr::read_lines(fastq_path)

######################################################## MAIN CODE
# get DNAStringSet of reads seqeuce
fastq_seq <- 
  fastq_file[seq(2, length(fastq_file), 4)] %>% 
  DNAStringSet(.)

# nucleotide frequency in the table
freq_tb <- 
  fastq_seq %>% 
  alphabetFrequency(., as.prob = T) %>% 
  .[, 1:4]

# get reads which have < 5% C or G nucleotides
converted_reads <- 
  fastq_file[seq(1, length(fastq_file), 4)][which(freq_tb[, 2] < 0.05 | freq_tb[, 3] < 0.05)] %>% 
  str_remove(., "^@")

# save to 
readr::write_lines(x = converted_reads, path = file.path(outpath, str_c(fastq_name, "converted_reads.2", "txt", sep = ".")), append = T)




