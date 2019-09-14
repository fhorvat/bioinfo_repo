#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: scales bigWig to RPMs
### DATE: Tue May 21 08:13:09 2019
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

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
bw_path <- args$raw_tracks
read_stats_path <- args$read_stats

######################################################## READ DATA
# read bigWig
bw <- rtracklayer::import(bw_path)

# read read counts
read_stats <- readr::read_delim(read_stats_path, delim = "\t")

######################################################## MAIN CODE
# get sample name
sample_id <- 
  basename(bw_path) %>% 
  str_remove(., "\\.bw$")

sample_id

# check if sample ID in read stats is equal to bigWig name
if(sample_id != read_stats$sample_id){
  stop("Something's not right - sample ID's in bigWig and stats table not matching. Please make sure your sample ID's match")
}

# scale bigWig
mcols(bw)$score <- (mcols(bw)$score) / (round((read_stats$mapped_minus_rDNA / 1e6), 6))

# write bw
rtracklayer::export(object = bw, 
                    con = str_c(sample_id, ".scaled.bw"))


