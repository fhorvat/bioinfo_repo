### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq")

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

# datasets path
datasets_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# list of datasets
datasets_list <- c("hamster_testis_Mov10l.8.5dpp.RNAseq/Analysis/expression.added_PIWIL3", 
                   "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq/Analysis/expression.added_PIWIL3.stranded", 
                   "hamster_testis_Mov10l.RNAseq/Analysis/expression.added_PIWIL3.stranded")

# list raw counts
raw_counts_list <- file.path(datasets_path,
                             datasets_list, 
                             "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.counts.txt")

# list FPKM tables
fpkm_tables_list <- file.path(datasets_path,
                              datasets_list, 
                              "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.FPKM.csv")

######################################################## READ DATA
# read all raw counts
raw_counts_tb <- purrr::map(raw_counts_list, function(path){

  # read 
  readr::read_delim(path, delim = "\t", comment = "#") %>%
    set_colnames(., basename(colnames(.))) %>% 
    dplyr::select(-c(Chr:Length)) %>% 
    dplyr::rename(gene_id = Geneid)
  
}) %>% 
  purrr::reduce(., left_join, by = "gene_id")

# read all FPKM tables
fpkm_tb <- purrr::map(fpkm_tables_list, function(path){
  
  # read 
  readr::read_csv(path)
  
}) %>% 
  purrr::reduce(., left_join) %>% 
  dplyr::select(gene_id, gene_name, coordinates, strand, gene_biotype, gene_description, everything())

######################################################## MAIN CODE
# write raw table
readr::write_csv(raw_counts_tb, 
                 file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.counts.all_samples.csv"))

# write FPKM table
readr::write_csv(fpkm_tb, 
                 file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.FPKM.all_samples.csv"))


