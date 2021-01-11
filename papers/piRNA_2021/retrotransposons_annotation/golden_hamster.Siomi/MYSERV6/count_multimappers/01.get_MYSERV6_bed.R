### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/MYSERV6/count_multimappers")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

######################################################## MAIN CODE
# subset MYSERV6-int
subset_gr <- 
  rmsk_clean %>% 
  dplyr::filter(repName == "MYSERV6-int") %>% 
  GRanges(.)

# save
rtracklayer::export.bed(subset_gr, file.path(outpath, "MYSERV6-int.rmsk.all.bed"))

