### INFO: 
### DATE: Thu Apr 23 22:11:27 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/miRNA_targets_expression/transcripts_expression")

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

library(tximport)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# transcripts info path
transcripts_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.transcriptInfo.csv$"), full.names = T)
  
# list salmon quantification paths
salmon_paths <- list.files(inpath, pattern = "quant.sf", full.names = T, recursive = T)
names(salmon_paths) <- 
  salmon_paths %>% 
  dirname(.) %>% 
  basename(.) %>% 
  str_remove(., "_quant")
  
######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read transcript-genes relation
transcripts_info <- 
  readr::read_csv(transcripts_info_path) %>% 
  dplyr::select(transcript_id, gene_id)

######################################################## MAIN CODE
# txi <- tximport(salmon_paths, type = "salmon", tx2gene = transcripts_info, ignoreTxVersion = T)

# read salmon quantifications 
salmon_list <- purrr::map(salmon_paths, function(path){
  
  # sample name
  sample_name <- 
    path %>% 
    dirname(.) %>% 
    basename(.) %>% 
    str_remove(., "_quant")
  
  # read and shape
  readr::read_delim(path, delim = "\t") %>% 
    dplyr::select(transcript_id = Name, TPM) %>% 
    dplyr::mutate(sample_id = sample_name)
  
}) %>% 
  bind_rows(.) 

# calculate mean TPM
tpm_tb <- 
  salmon_list %>% 
  dplyr::group_by(transcript_id) %>% 
  dplyr::summarise(TPM = mean(TPM)) %>% 
  dplyr::mutate(transcript_id = str_remove(transcript_id, "\\..*"))

# for each gene get transcript with higher expression
gene_tb <- 
  transcripts_info %>% 
  left_join(., tpm_tb, by = "transcript_id") %>% 
  dplyr::select(gene_id, transcript_id, TPM) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(TPM_rank = rank(-TPM, ties.method = "random")) %>%
  dplyr::filter(TPM_rank == 1) %>% 
  dplyr::select(-TPM_rank)
  
# save 
gene_tb %>% 
  readr::write_csv(., file.path(outpath, "ensembl.93.Sscrofa11.1.20180920.UCSCseqnames.gene_transcript_expression.pig_GV.csv"))


