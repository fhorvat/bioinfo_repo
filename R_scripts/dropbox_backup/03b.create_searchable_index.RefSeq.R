### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/RefSeq")

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
fasta_path <- list.files(inpath, ".*fna.gz", full.names = T)

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
  tibble(gene_id = str_remove(fasta_tb, ">") %>% str_remove(., " .*"),
         gene_symbol = str_remove(fasta_tb, ".*\\(") %>% str_remove(., "\\).*")) %>% 
  dplyr::mutate(gene_id_2 = gene_id)

# save index
geneID_tb_out <-
  geneID_tb %>% 
  dplyr::filter(!is.na(gene_id), !is.na(gene_symbol)) %T>% 
  readr::write_delim(., path = file.path(outpath, "GCF_000349665.1_MesAur1.0_rna.geneID_to_geneSymbol.txt"), delim = " ", col_names = F)
