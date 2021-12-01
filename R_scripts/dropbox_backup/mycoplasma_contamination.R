### INFO: 
### DATE: Thu Sep 20 23:04:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/mycoplasma_contamination")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# logs paths
"/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/mycoplasma_contamination"
"/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Zuzka_3T3_PAPD7/Data/Raw/mycoplasma_contamination"
"/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Raw/mycoplasma_contamination"
"/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Raw/mycoplasma_contamination"
"/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/mouse/Data/Raw/mycoplasma_contamination"
"/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Raw/mycoplasma_contamination"

######################################################## READ DATA

######################################################## MAIN CODE
### grep logs from dir
# smallRNA seq
smallRNA <- 
  bind_rows(data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Zuzka_3T3_PAPD7/Data/Raw/mycoplasma_contamination/*k27.log', 
                            header = F), 
          data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Raw/mycoplasma_contamination/*k27.log', 
                            header = F)) %>% 
  tibble::as.tibble(.) %>%
  dplyr::mutate(V1 = basename(V1) %>% str_remove(., "\\.k27.log.*$")) %>%
  dplyr::select(sample_id = V1, reads = V2, bases = V3) %T>%
  readr::write_csv(., file.path(outpath, "mycoplasma_contamination.smallRNAseq.csv"))
  
# mESCs long RNAseq
longRNA <- 
  bind_rows(data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Raw/mycoplasma_contamination/*k27.log', 
                            header = F), 
          data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/mouse/Data/Raw/mycoplasma_contamination/*k27.log', 
                            header = F)) %>% 
  tibble::as.tibble(.) %>%
  dplyr::mutate(V1 = basename(V1) %>% str_remove(., "\\.k27.log.*$")) %>%
  dplyr::select(sample_id = V1, reads = V2, bases = V3) %T>%
  readr::write_csv(., file.path(outpath, "mycoplasma_contamination.longRNAseq.csv"))

# lncKO RNAseq
bind_rows(data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Raw/mycoplasma_contamination/*k27.log', 
                            header = F)) %>% 
  tibble::as.tibble(.) %>%
  dplyr::mutate(V1 = basename(V1) %>% str_remove(., "\\.k27.log.*$")) %>%
  dplyr::select(sample_id = V1, reads = V2, bases = V3) %T>%
  readr::write_csv(., file.path(outpath, "mycoplasma_contamination.lncKO.RNAseq.csv"))

# positive and negative control
bind_rows(data.table::fread('grep "Contaminants:" /common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/mycoplasma_contamination/*k27.log', 
                            header = F)) %>% 
  tibble::as.tibble(.) %>%
  dplyr::select(sample_id = V1, reads = V2, bases = V3) %>% 
  dplyr::mutate(sample_id = basename(sample_id) %>% str_remove(., "\\.k27.log.*$"), 
                sample_id = ifelse(sample_id == "SRR944282", "positive_control", "negative_control")) %>% 
  asdf
