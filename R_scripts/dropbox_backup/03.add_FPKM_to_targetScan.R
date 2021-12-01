### INFO: 
### DATE: Thu Apr 25 13:27:01 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/last_exon.let7_counts")

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

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
gene_info_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*geneInfo.csv", full.names = T)

# transcript info path
transcript_info_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.transcriptInfo.csv")

# oocyte expressed genes
oocyte_fpkm_path <- file.path(inpath, "ensembl.93.mouse_B6.mm10.avgFPKM.RDS")

# target scan path
targetScan_path <- file.path(inpath, "TargetScan7.1__let-7-5p_miR-98-5p.predicted_targets.txt")

######################################################## READ DATA
# read gene info
genes_info <- readr::read_csv(gene_info_path)

# read transcripts info
transcripts_info <- readr::read_csv(transcript_info_path)

# read oocyte expressed genes
oocyte_fpkm <-readRDS(oocyte_fpkm_path)

# read targetScan
targetScan_tb <- readr::read_delim(targetScan_path, delim = "\t")

######################################################## MAIN CODE
# clean FPKM table
oocyte_fpkm_tidy <- 
  oocyte_fpkm %>% 
  as.tibble(.) %>%
  dplyr::select(gene_id, GV_WT)

# add FPKMs to transcipts IDs
transcripts_info_fpkm <- 
  transcripts_info %>% 
  left_join(., oocyte_fpkm_tidy, by = "gene_id")

# add FPKM to targetScan table
targetScan_fpkm <- 
  targetScan_tb %>% 
  dplyr::mutate(transcript_id = str_remove(`Representative transcript`, "\\.[0-9]+$")) %>% 
  left_join(., transcripts_info_fpkm, by = "transcript_id") %>% 
  dplyr::select(`Target gene`, CNOT6L.GV_WT.FPKM = GV_WT, `Representative transcript`, gene_id.ens93 = gene_id, everything()) %>% 
  dplyr::select(-c(transcript_id, gene_name)) %>% 
  dplyr::filter(CNOT6L.GV_WT.FPKM >= 2) %T>%
  write_csv(., file.path(outpath, basename(targetScan_path) %>% str_replace(., "\\.txt", ".GV_FPKM.2_FPKM_cut.20190425.csv")))
  

