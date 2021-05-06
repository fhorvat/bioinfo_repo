### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/Liftoff/BCM_Maur_2.0.GCA_017639785.1/RefSeq/hub")

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
gtf_path <- file.path(inpath, "..")
gtf_path <- file.path(gtf_path, "GCF_017639785.1_BCM_Maur_2.0_genomic.Siomi.liftoff.fixed.gff")

######################################################## READ DATA
# read gtf using rtracklayer
gtf_gr <- rtracklayer::import.gff(gtf_path)

######################################################## MAIN CODE
# save as gff3
rtracklayer::export.gff3(gtf_gr, file.path(outpath, "GCF_017639785.1_BCM_Maur_2.0_genomic.Siomi.liftoff.gff3"))
