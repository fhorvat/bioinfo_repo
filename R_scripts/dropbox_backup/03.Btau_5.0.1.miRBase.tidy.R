### INFO: 
### DATE: Mon Feb 03 07:54:27 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/cow/Btau_5.0.1.GCA_000003205.6")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# fasta index path
fai_path <- file.path(inpath, "Btau_5.0.1.fa.fai")

# assembly report path
assembly_path <- file.path(inpath, "GCA_000003205.6_Btau_5.0.1_assembly_report.txt")

# miRBase gff3
mirbase_path <- file.path(inpath, "miRBase.22.Btau_5.0.1.20200202.gff3")

######################################################## READ DATA
# read fasta index
fai_tb <- 
  readr::read_delim(file = fai_path, delim = "\t", col_name = c("genBank", "width", "tmp1", "tmp2", "tmp3")) %>%
  dplyr::select(genBank)

# read assembly report
assembly_tb <- readr::read_delim(file = assembly_path, delim = "\t", comment = "#", 
                                 col_names = c("sequence_name", "sequence_role", "assigned_molecule",
                                               "assigned_molecule_location_type", "genBank_accn",
                                               "relationship", "refSeq_accn", "assembly_unit",
                                               "sequence_length", "UCSC_style_name"))

# read miRBase gff3
mirbase_tb <- 
  rtracklayer::import(mirbase_path) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(seqnames = as.character(seqnames))

######################################################## MAIN CODE
# clean assembly report
assembly_tb_tidy <- 
  assembly_tb %>% 
  dplyr::select(seqnames = sequence_name, genBank = genBank_accn) %>% 
  dplyr::mutate(seqnames = str_replace(seqnames, "Chromosome", "chr"))

# change chromosome levels in miRBase .gff
mirbase_tb_tidy <- 
  mirbase_tb %>% 
  dplyr::left_join(., assembly_tb_tidy, by = "seqnames") %>% 
  dplyr::select(-seqnames) %>% 
  dplyr::select(seqnames = genBank, everything())

# save as .gff3
mirbase_tb_tidy %>% 
  GRanges(.) %T>%
  rtracklayer::export.gff3(., con = file.path(outpath, "miRBase.22.Btau_5.0.1.20200202.genBankseqnames.gff3"))