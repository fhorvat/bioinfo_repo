### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/Papd7_KO/polyA_tails")

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

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- file.path(genome_dir, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.gtf.gz")

# original gene info path
geneInfo_path <- file.path(genome_dir, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read geneInfo
geneInfo <- readr::read_csv(geneInfo_path)

######################################################## MAIN CODE
# get 3p UTR, resize to last 10 bases
gtf_3UTR <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "three_prime_utr") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_name, gene_biotype, transcript_id) %>% 
  GRanges(.) %>% 
  GenomicRanges::resize(., 10, fix = "end")

# save as table
readr::write_csv(as_tibble(gtf_3UTR), file.path(outpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.3pUTR.annotation.csv"))

# save as .bed
gtf_3UTR %T>% 
  rtracklayer::export.bed(., file.path(outpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.3pUTR.annotation.bed"))

