### INFO: 
### DATE: Mon Jul 08 12:52:07 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/UCSC_maf")

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
library(Biostrings)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# strains path
strains_path <- file.path(inpath, "mouse_strains")

# multiz60way rat
rat_path <- file.path(inpath, "rat_multiz60way")

# lnc1 mouse strains bed path, get names
lnc1_bed_paths <- list.files(c(strains_path, rat_path), pattern = ".*\\.mm10\\.lnc1_3prime\\.bed", full.names = T)

######################################################## READ DATA
# read, reduce and save bed
purrr::map(lnc1_bed_paths, function(path){
    
    # filter path, import and get ranges of bed
    maf_bed <- 
      rtracklayer::import.bed(path) %>% 
      reduce(.) %>% 
      range(.)
    
    # .maf is 1-based, start should be shifted by one base
    start(maf_bed) <- start(maf_bed) - 1
    
    # save as .bed
    rtracklayer::export.bed(maf_bed, file.path(outpath, str_replace(basename(path), "\\.mm10\\.lnc1_3prime\\.bed$", ".mm10.lnc1_3prime.top_result.bed")))
    
  })

######################################################## MAIN CODE


