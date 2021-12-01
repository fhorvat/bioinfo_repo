### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/blast_ORFs_L1")

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

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(inpath, "../expression")

# table with annotation and coordinates
results_tb_path <- file.path(inpath, "L1_full_length_manual_200707.consensus.ORFs_blast_hits.csv")

######################################################## READ DATA
# read table with annotation and coordinates of results
results_tb <- readr::read_csv(results_tb_path) 

######################################################## MAIN CODE
# create out name
out_name <- 
  results_tb_path %>% 
  basename(.) %>% 
  stringr::str_replace(., "\\.csv$", ".saf")

# check if the file name is the same as the input table, stop if it is
if(file.path(outpath, out_name) != results_tb_path){
  
  # save as SAF
  results_saf <- 
    results_tb %>% 
    tidyr::separate(hit_coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
    dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
    readr::write_delim(., file.path(outpath, out_name), delim = "\t")
  
}


