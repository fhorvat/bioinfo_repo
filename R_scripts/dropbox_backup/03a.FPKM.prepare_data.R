### INFO: 
### DATE: Tue Oct 29 22:03:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test")

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

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 89

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# RDS path
rds_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read assay from .RDS
se <- readRDS(rds_path)
  
######################################################## MAIN CODE
# get only genes on chromosome 5
chr5_genes <- 
  genes_info %>% 
  dplyr::filter(seqnames == "chr5") %$%
  gene_id

# get total length of all exons for each transcript
gene_info <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), total_exon_length_sum = .) %>% 
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::filter(seqnames %in% c("chr5", "chr6")) %>% 
  dplyr::select(-c(seqnames:strand)) %T>%
  write_csv(., "ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.chr5_chr6.csv")

# filter se
se <- se[rownames(se) %in% chr5_genes, ]
assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %T>%
  write_csv(., "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.chr5.csv")

