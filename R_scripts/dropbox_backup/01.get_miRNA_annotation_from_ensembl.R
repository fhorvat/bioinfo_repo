### INFO: gets miRNA annotation
### DATE: Fri Jun 04 16:23:06 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/annotation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

######################################################## READ DATA
### get miRNA annotations
# set genome list
# genome_list <- c("cow.bosTau9", "ghamster.mesAur1", "human.hg38", "mouse.mm10", "pig.susScr11")
genome <- "ghamster.mesAur1"

## paths
# genome path
genome_dir <- file.path(inpath, genome)

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)


## read files
# read gene info
genes_info <- readr::read_csv(genes_info_path)

# read exons
exons_gr <-
  readRDS(file = exons_path) %>%
  tibble::as_tibble(.) %>%
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name, gene_biotype, gene_description), by = "gene_id") %>%
  dplyr::filter(gene_biotype == "miRNA") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_name, gene_biotype, gene_description) %>%
  dplyr::mutate(gene_biotype = replace(gene_biotype, gene_biotype == "miRNA", "miRNA.ENSEMBL"), 
                type = "miRNA") %>%
  dplyr::rename(ID = gene_id, Name = gene_name) %>%  
  dplyr::select(-c(gene_biotype, gene_description)) %>% 
  GenomicRanges::GRanges(.)

# save as gff
rtracklayer::export.gff3(exons_gr, con = file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.miRNAs.gff3"))
