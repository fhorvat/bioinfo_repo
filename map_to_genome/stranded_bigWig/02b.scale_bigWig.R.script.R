#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: scales bigWig to RPMs
### DATE: Tue May 21 08:13:09 2019
### AUTHOR: Filip Horvat

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

library(rtracklayer)

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

# list of read stats, genome logs and merged logs
bw_path <- args$bw_path
scale_factor <- as.numeric(args$scale_factor)

######################################################## READ DATA
# read bigWig
bw <- rtracklayer::import(bw_path)

######################################################## MAIN CODE
# get sample name
sample_name <-
  basename(bw_path) %>%
  str_remove(., "\\.bw$")

# scale bigWig
mcols(bw)$score <- mcols(bw)$score / scale_factor

# write bw
rtracklayer::export(object = bw,
                    con = str_c(sample_name, ".scaled.bw"))
