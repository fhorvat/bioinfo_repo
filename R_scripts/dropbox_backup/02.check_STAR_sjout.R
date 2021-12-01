### INFO: read SJ,out.tab from CNOT6L GV WT data and look for the splice junction 
### DATE: Thu Mar 08 21:39:48 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10_noMultimapFilter")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "read_SJout.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read SJ.out from CNOT6L GV WT data
sjout <- 
  list.files(path = outpath, pattern = "s_GV_WT_r1.SJ.out.tab") %>% 
  read_SJout(.) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
# find overlaps with splice acceptor site chr15:73134849-73134852
splice_acceptor <- GenomicRanges::GRanges(seqnames = "chr15", ranges = IRanges::IRanges(start = 73134849, end = 73134852), strand = "*")
splice_overlaps <- findOverlaps(splice_acceptor, sjout, ignore.strand = T)
sjout[subjectHits(splice_overlaps)]
