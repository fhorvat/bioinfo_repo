### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Siomi_assembly")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# liftoff gtf path
gtf_path <- file.path(inpath, "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff")

######################################################## READ DATA
# read gtf using rtracklayer
gtf_gr <- rtracklayer::import.gff(gtf_path)

# # read gtf as table
# gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# change seqnames
seqlevels(gtf_gr) <- str_c(seqlevels(gtf_gr), "|arrow")

# save as gff3
rtracklayer::export.gff3(gtf_gr, file.path(outpath, "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff3"))

######################################################## MAIN CODE
