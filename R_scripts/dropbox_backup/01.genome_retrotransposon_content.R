### INFO: 
### DATE: Wed Aug 18 19:08:48 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference"
mouse_path <- file.path(genome_dir, "mouse/mm10.GRCm38.GCA_000001635.2")
human_path <- file.path(genome_dir, "human/hg38.GRCh38.GCA_000001405.15")

# get repeatMasker path
rmsk_path_list <- list.files(c(mouse_path, human_path), "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

# get fasta index paths
fai_path_list <- list.files(c(mouse_path, human_path), ".*\\.fa\\.fai", full.names = T)

######################################################## READ DATA
# read repeatMasker tables
rmsk_list <- purrr::map(rmsk_path_list, function(path){
  
  # read
  readr::read_delim(path, delim = "\t")
  
}) %>%
  set_names(str_extract(rmsk_path_list, "human|mouse"))

# read fasta index files
fai_list <- purrr::map(fai_path_list, function(path){
    
    # read
    readr::read_delim(path, delim = "\t", col_names = c("seqnames", "length", "tmp1", "tmp2", "tmp3"))
    
  }) %>%
  set_names(str_extract(fai_path_list, "human|mouse"))

######################################################## MAIN CODE
# get total genome lengths
genome_length <- purrr::map(c("human", "mouse"), function(animal){
  
  # get the table, sum the lengths
  genome_tb <- 
    fai_list[[animal]] %>% 
    dplyr::summarise(genome_length = sum(length)) %>% 
    dplyr::mutate(animal = animal)
  
  # return
  return(genome_tb)
  
}) %>% 
  dplyr::bind_rows(.)
  

# get lengths of retrotransposons
retro_length <- purrr::map(c("human", "mouse"), function(animal){
  
  # get repeatMasker table, filter
  rmsk_tb <- 
    rmsk_list[[animal]] %>% 
    dplyr::filter(repClass %in% c("LINE", "SINE", "LTR", "SINE?", "LINE?", "LTR?"))
  
  # create GRanges, reduce iregardless of strand
  rmsk_reduced <- 
    rmsk_tb %>% 
    GRanges(.) %>% 
    reduce(., ignore.strand = T)
  
  # get the total length of all retrotranspsons
  retro_tb <- 
    rmsk_reduced %>% 
    tibble::as_tibble(.) %>% 
    dplyr::summarise(retro_length = sum(width)) %>% 
    dplyr::mutate(animal = animal)
  
  # return
  return(retro_tb)
  
}) %>% 
  dplyr::bind_rows(.)


# join tables, get the percentage of retrotransposon derived sequences
genome_retro_tb <- 
  genome_lengths %>% 
  dplyr::left_join(., retro_length, by = "animal") %>% 
  dplyr::mutate(retro_perc = round((retro_length / genome_length), 3) * 100)
  