### INFO: 
### DATE: Fri Jun 04 18:41:42 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2/miRBase_liftoff")

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
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/DB/genome_reference/cow/Btau_5.0.1.GCA_000003205.6"

# set outpath
outpath <- getwd()

# miRBase gff3
mirbase_path <- file.path(inpath, "miRBase.22.Btau_5.0.1.20200202.genBankseqnames.gff3")
mirbase_out_path <- file.path(outpath, "miRBase.22.Btau_5.0.1.20200202.genBankseqnames.gff3")

# check 
if(mirbase_path == mirbase_out_path){
  stop("Carful! You'll overwrite your original .gff3 file")
}

######################################################## READ DATA
# read miRBase gff3
mirbase_tb <-
  rtracklayer::import(mirbase_path) %>%
  as_tibble(.) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) 

######################################################## MAIN CODE
# rename
mirbase_tb %>% 
  dplyr::rename(Parent = Derives_from) %>% 
  GRanges(.) %T>%
  rtracklayer::export.gff3(., con = file.path(outpath, "miRBase.22.Btau_5.0.1.20200202.genBankseqnames.gff3"))

