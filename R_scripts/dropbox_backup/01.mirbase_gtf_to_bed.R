### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May/Analysis/Dicer_KO_reads")

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
gtf_path <- list.files(path = genome_dir, pattern = "miRBase.*\\.gff3", full.names = T)

######################################################## READ DATA
# read miRBase gtf
mirna_gr <- rtracklayer::import.gff(con = gtf_path)

######################################################## MAIN CODE
# read and clean miRBase miRNA annoation
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("ID")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

# save as bed
rtracklayer::export.bed(mirna_gr, file.path(outpath, "miRBase.22.mm10.20181605.mature_miRNA.bed"))
