### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/ENSEMBL")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# fasta path
fasta_path <- list.files(inpath, ".*fa.gz", full.names = T)

######################################################## READ DATA
# read fasta headers
fasta_head <- readr::read_lines(fasta_path)

######################################################## MAIN CODE
# filter
fasta_tb <- 
  fasta_head %>% 
  .[str_detect(., ">")]

# create table
geneID_tb <- 
  fasta_tb %>% 
  stringr::str_split(., pattern = " ") %>% 
  purrr::map(., function(x) tibble(gene_id = x[str_detect(x, "gene:")], 
                                   gene_symbol = x[str_detect(x, "gene_symbol:")])) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "gene:"), 
                gene_symbol = str_remove(gene_symbol, "gene_symbol:"), 
                gene_id_2 = gene_id)

# save index
geneID_tb_out <-
  geneID_tb %>% 
  dplyr::filter(!is.na(gene_id), !is.na(gene_symbol)) %T>% 
  readr::write_delim(., path = file.path(outpath, "ensembl.99.MesAur1.0.20200609.cdna.all.geneID_to_geneSymbol.txt"), delim = " ", col_names = F)
