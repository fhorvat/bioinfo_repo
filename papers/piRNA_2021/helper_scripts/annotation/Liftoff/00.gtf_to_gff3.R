### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/Liftoff/MesAur1/ENSEMBL")

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
gtf_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"
gtf_path <- file.path(gtf_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")

######################################################## READ DATA
# read gtf using rtracklayer
gtf_gr <- rtracklayer::import.gff(gtf_path)

# save as gff3
rtracklayer::export.gff3(gtf_gr, file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gff3"))
