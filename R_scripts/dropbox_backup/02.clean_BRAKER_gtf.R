### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hubs/golden_hamster.Siomi/files/denovo_annotation")

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

# braker gtf path
braker_gtf_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/BRAKER/braker"
braker_gtf_path <- file.path(braker_gtf_path, "braker.gtf")

######################################################## READ DATA
# read gtf using rtracklayer
braker_gtf <- rtracklayer::import.gff(braker_gtf_path)

# read gtf as table
braker_gtf_tb <- read_delim(file = braker_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# write.table(x = ., file = file.path(outpath, gtf_name), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

######################################################## MAIN CODE
# add gene_id to genes
braker_gtf_tb_clean <- 
  braker_gtf_tb %>% 
  dplyr::mutate(X9 = ifelse(X3 == "gene", str_c("gene_id ", "\"", X9, "\""), X9))

# save
write.table(x = braker_gtf_tb_clean, file = file.path(outpath, "braker.clean.gtf"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


