### INFO: gets 3'UTR sequences, longest for each gene
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/miR-205_pig/oocyte_seed_expression")

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

library(BSgenome.Sscrofa.UCSC.susScr11)
library(Biostrings)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# set sequence pattern
seed_pattern_list <- 
  c("ATGAAGG", "AATGAAGG") %>% 
  set_names(., c("7mer", "8mer"))
  
# loop for 7mer and 8mer
purrr::map(names(seed_pattern_list), function(name){
  
  # match pattern to the genome
  seed_match <- Biostrings::vmatchPattern(pattern = seed_pattern_list[[name]], subject = BSgenome.Sscrofa.UCSC.susScr11)
  
  # set names
  names(seed_match) <- str_c(seqnames(seed_match), ":", start(seed_match), "-", end(seed_match))
  
  # export as .gtf
  export.gff3(object = seed_match, 
              con = file.path(outpath, 
                              str_c("miR-205", "susScr11", name, "seed_match", "gff", sep = ".")))
  
  
  # return 
  return(name)
  
})


