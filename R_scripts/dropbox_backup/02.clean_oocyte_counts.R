### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Documentation")

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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"
experiment_paths <- c("hamster_oocyte_Mov10l.RNAseq/Analysis/expression.added_PIWIL3")
experiment_paths <- file.path(base_path, experiment_paths, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.counts.txt")

######################################################## READ DATA
# read count tables
counts_tb <- purrr::map(experiment_paths, function(path){
  
  readr::read_delim(path, delim = "\t", comment = "#") %>%
    set_colnames(., basename(colnames(.))) %>% 
    dplyr::select(-c(Chr:Length)) %>%
    dplyr::rename(gene_id = Geneid)
  
  
}) %>% 
  purrr::reduce(., left_join, by = "gene_id") %>% 
  magrittr::set_colnames(., str_remove(colnames(.), "\\.bam$"))

######################################################## MAIN CODE
# tidy and save table
counts_tidy <- 
  counts_tb %>% 
  tidyr::pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "count") %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT|HET|KO"), 
                replicate = str_extract(sample_id, "(?<=r)[1-9]+(?=\\.)") %>% as.numeric(.)) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "HET", "KO"))) %>% 
  dplyr::arrange(genotype, replicate) %>% 
  tidyr::pivot_wider(id_cols = gene_id, names_from = sample_id, values_from = count) %T>%
  readr::write_delim(., file = file.path(outpath, "ensembl.99.MesAur1.0.Piwil3_Tex101.from_RefSeq.raw_counts.oocyte.txt"))



