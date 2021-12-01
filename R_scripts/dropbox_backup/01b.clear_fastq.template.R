### INFO: read fastq, discard reads which have "unfixable" in ID, strip "cor" from ID, write
### DATE: 7. 12. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(ShortRead)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# get one fastq file
fastq_file <- "%FILE"

######################################################## MAIN CODE
# get fastq name
fastq_name <- 
  basename(fastq_file) %>% 
  str_replace(fastq_file, "cor", "cor2")

# set filter - "unfixable" in ID
filt <- idFilter(regex = "unfixable", fixed = FALSE, exclude = T)

# filter fastq
fastq_file <- fastq_file[filt(fastq_file)]

# remove "cor" from read ID
fastq_filt@id[] <- 
  str_replace(ShortRead::id(fastq_filt), " cor", "") %>% 
  BStringSet(.)

# write corrected .fastq
ShortRead::writeFastq(object = fastq_filt, file = fastq_name, compress = F)
