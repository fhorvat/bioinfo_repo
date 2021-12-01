### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/3pUTR_hexamers")

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
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gtf path
gtf_path <- file.path(genome_dir, "ensembl.100.Sscrofa11.1.20200616.UCSCseqnames.gtf.gz")

# original gene info path
geneInfo_path <- file.path(genome_dir, "ensembl.100.Sscrofa11.1.20200616.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read geneInfo
geneInfo <- readr::read_csv(geneInfo_path)

######################################################## MAIN CODE
# get 3p UTR
gtf_3UTR <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "three_prime_utr") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_name, gene_biotype, transcript_id) %>% 
  GRanges(.)

# save as table
readr::write_csv(as_tibble(gtf_3UTR), file.path(outpath, "pig.susScr11.ensembl.100.Sscrofa11.1.20200616.UCSCseqnames.3pUTR.annotation.csv"))

# save as .bed
gtf_3UTR %T>% 
  rtracklayer::export.bed(., file.path(outpath, "pig.susScr11.ensembl.100.Sscrofa11.1.20200616.UCSCseqnames.3pUTR.annotation.bed"))

