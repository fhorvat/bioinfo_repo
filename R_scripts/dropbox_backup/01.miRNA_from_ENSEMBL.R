### INFO: 
### DATE: Sun Feb 02 09:51:30 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("")

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

# set ensembl version
ensembl_version <- 96

# genome path
genome_dir <- "/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
# miRNA 
mirna_ensembl <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "miRNA")

# save as gff3
rtracklayer::export.gff3(mirna_ensembl, file.path(genome_dir, basename(genes_info_path) %>% str_replace(., "geneInfo.csv$", "miRNA.gff3")))
