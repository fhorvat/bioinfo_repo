### INFO: get connection between ensembl and refSeq scaffolds in GCA_000349665.1 MesAur1.0 assembly, rename seqnames in ensembl gtf
### DATE: Sat Mar 10 21:20:05 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/reference/hamster_golden/mesAur1/vfranke")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(biomaRt)

######################################################## PATH VARIABLES
outpath <- getwd()

# scaffold names path
scaffold_path <- "/common/WORK/fhorvat/reference/hamster_golden/mesAur1/vfranke/scaffold_names.txt"
  
# ensembl gtf path
ensembl_gtf <- "/common/WORK/fhorvat/reference/hamster_golden/mesAur1/ensembl/ensembl.MesAur1.0.91.20180309.gtf.gz"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "GffToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read table with scaffold names
scaffold_df <- 
  read_delim(file = scaffold_path, delim = "\t") %>% 
  magrittr::set_colnames(., str_replace_all(colnames(.), " |#", ""))

# read ensembl gtf
ensembl_df <- read_delim(file = ensembl_gtf, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# add ensembl annotation, change mitochondrial scaffold, save
scaffold_df %<>% 
  dplyr::mutate(ensembl = GenBankAccession.version) %>% 
  dplyr::mutate(ensembl = replace(ensembl, ensembl == "EU660218.1", "MT")) %>% 
  dplyr::select(RefSeq_accession = RefSeqAccession.version, ensembl_name = ensembl) %T>% 
  readr::write_delim(., path = file.path(outpath, "GCF_000349665.1_MesAur1.0.ensembl_refSeq.annotation.txt"), delim = "\t")

# change seqnames in ensembl .gtf from GenBank to refSeq
ensembl_df_rename <- 
  ensembl_df %>% 
  dplyr::left_join(., scaffold_df, by = c("X1" = "ensembl_name")) %>% 
  dplyr::mutate(X1 = RefSeq_accession) %>% 
  dplyr::select(X1, everything(), -RefSeq_accession) %>% 
  readr::write_delim(., path = file.path(outpath, "ensembl.refSeq_annotation.MesAur1.0.91.20180309.gtf.gz"), delim = "\t", col_names = F)
