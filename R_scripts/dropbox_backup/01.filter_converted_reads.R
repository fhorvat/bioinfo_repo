### INFO: 
### DATE: Wed Sep 30 13:28:24 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/bisulfite_sequencing/test_run/conversion_status")

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

library(ShortRead)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# data path
data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Raw/Links"

# list fastq files
fastq_list <- list.files(data_path, "s_mouse.*\\.txt\\.gz", full.names = T)
  
######################################################## READ DATA
### sample 1 million reads in first and second read in a pair
# sample fastq
fastq_sample_1 <- FastqSampler(fastq_list[1])

# get reads
fastq_reads_1 <- yield(fastq_sample_1)

# get sequences
fastq_seq_1 <- sread(fastq_reads_1)
fastq_id_1 <- id(fastq_reads_1)


### second read in a pair
# sample fastq
fastq_sample_2 <- FastqSampler(fastq_list[2])

# get reads
fastq_reads_2 <- yield(fastq_sample_2)

# get sequences
fastq_seq_2 <- sread(fastq_reads_2)
fastq_id_2 <- id(fastq_reads_2)

######################################################## MAIN CODE
# in the table
freq_tb <- 
  fastq_seq_2 %>% 
  alphabetFrequency(., as.prob = T) %>% 
  .[, 1:4] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(read_id = as.character(fastq_id)) %>% 
  dplyr::select(read_id, everything())

# filter
freq_tb %>% 
  dplyr::filter(C < 0.05) %>% 
  nrow(.) %>% 
  magrittr::divide_by(., length(fastq_seq)) %>% 
  magrittr::multiply_by(., 100)

freq_tb %>% 
  dplyr::filter(G < 0.05) %>% 
  nrow(.) %>% 
  magrittr::divide_by(., length(fastq_seq)) %>% 
  magrittr::multiply_by(., 100)

