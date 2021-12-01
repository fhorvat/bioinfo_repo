### INFO: 
### DATE: Sat Oct 05 13:47:57 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### IN AND OUT
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### COMMAND LINE ARGUMENTS
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

experiment <- args$experiment
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)
dataset_path <- args$dataset_path
feature_coordinates <- args$feature_coordinates


### DATASET
# bigWig path
coverage_path <- 
  list.files(dataset_path, pattern = ".*\\.bw", full.names = T) %>% 
  .[!str_detect(., "scaled")]


######################################################## READ DATA
# read LINE1 coordinates
line1_gr <- readRDS(feature_coordinates)
line1_tb <- 
  as_tibble(line1_gr) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

# read coverage from bigWig
coverage_list <- purrr::map(coverage_path, function(path){
  
  # read bigWig file
  coverage <- 
    rtracklayer::import(path) %>% 
    .[mcols(.)$score > 0] %>% 
    reduce(., ignore.strand = T)
  
}) %>% 
  set_names(., coverage_path %>% basename(.) %>% str_remove_all(., "\\.bw$"))

######################################################## MAIN CODE
### find coverage for each LINE1 element
line1_coverage_list <- purrr::map(names(coverage_list), function(sample_id){
  
  # get individual sample coverage
  coverage <- coverage_list[[sample_id]]
  
  # find overlaps between two GRanges
  hits <- findOverlaps(line1_gr, coverage, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  line1_coverage <- extractList(coverage, as(hits, "List"))
  
  # intersect with LINE1 coordinates
  line1_coverage <- pintersect(line1_coverage, line1_gr, ignore.strand = T)
  
  # set names
  names(line1_coverage) <- mcols(line1_gr)$rmsk_id 
  
  
  ### find longest stretch and percentage of coverage for each LINE1
  # unlist
  line1_coverage_tb <- 
    line1_coverage %>% 
    unlist(.)
  mcols(line1_coverage_tb)$rmsk_id <- names(line1_coverage_tb)
  
  # convert to data.table, summarize 
  line1_coverage_sum <- 
    line1_coverage_tb %>% 
    as_tibble(.) %>% 
    dplyr::group_by(rmsk_id) %>% 
    dplyr::summarise(coverage_total = sum(width), 
                     coverage_longest = max(width)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., line1_tb %>% dplyr::select(rmsk_id, width), by = "rmsk_id") %>% 
    dplyr::mutate(coverage_ratio = round((coverage_total / width), 3)) %>% 
    dplyr::select(rmsk_id, coverage_ratio, longest_stretch = coverage_longest) %>% 
    set_colnames(., c("rmsk_id", str_c(c("coverage_ratio", "longest_stretch"), sample_id, sep = ".")))
  
  # return 
  return(line1_coverage_sum)
  
})

# join to one table
line1_coverage_tb <- 
  purrr::reduce(c(list(line1_tb), line1_coverage_list), left_join, by = "rmsk_id") %>% 
  dplyr::mutate_at(vars(starts_with("coverage_ratio")), list(~ifelse(is.na(.), 0, .))) %>% 
  dplyr::mutate_at(vars(starts_with("longest_")), list(~ifelse(is.na(.), 0, .))) %T>%
  readr::write_csv(., path = file.path(outpath, str_c((feature_coordinates %>% basename(.) %>% str_remove(., "\\.GRanges\\.RDS$")), 
                                                      experiment, 
                                                      "coverage.csv", 
                                                      sep = ".")))
