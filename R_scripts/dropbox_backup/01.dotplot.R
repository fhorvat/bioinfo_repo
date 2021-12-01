### INFO: 
### DATE: Fri May 11 16:13:56 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/sequences/pairwise_alignments")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(dotplot)
library(Biostrings)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# sequences path
seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/sequences/mammals_plusStrand"

######################################################## READ DATA
# read sequences
sequences <- Biostrings::readDNAStringSet(filepath = list.files(path = seq_path, pattern = "*.fasta", full.names = T))

######################################################## MAIN CODE
# get intron3_4
intron_filt <- "intron3_4"
sequences_filt <- sequences[str_detect(names(sequences), intron_filt)]

# get characters
sequences_string <- 
  as.character(sequences_filt) %>%
  magrittr::set_names(., stringr::str_remove(names(.), " .*$"))
  
# set species to compare with mouse, get mous
species <- "rnorvegicus"
species_seq <- sequences_string[str_detect(names(sequences_string), species)]
mmusculus_seq <- sequences_string[str_detect(names(sequences_string), "mmusculus")]
  
# do dotplot
dotPlotg(species_seq, mmusculus_seq, wsize = 7) + 
  theme_bw() + 
  labs(x = "mmusculus", y = species) +
  ggsave(filename = file.path(outpath, str_c("dotplot.", "mmusculus_vs_", species, ".intron3_4.png")), 
         height = 10, width = 10)

