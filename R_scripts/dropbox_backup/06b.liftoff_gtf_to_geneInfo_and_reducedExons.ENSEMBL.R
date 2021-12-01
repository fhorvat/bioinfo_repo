### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# gtf path
gtf_path <- file.path(genome_dir, "annotation/Liftoff/MesAur1/ENSEMBL", "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.gff")

# original gene info path
geneInfo_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"
geneInfo_path <- file.path(geneInfo_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read geneInfo
geneInfo <- readr::read_csv(geneInfo_path)

######################################################## MAIN CODE
# remove exons without fully mapped genes
filt_list <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %$% 
  gene_id %>% 
  unique(.)

# filter gtf
gtf_gr <- gtf_gr[(mcols(gtf_gr)$gene_id %in% filt_list)]

# remove low identity and partially mapped genes
filt_list <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(!is.na(low_identity) | !is.na(partial_mapping)) %$%
  gene_id %>% 
  unique(.)

# filter gtf
gtf_gr <- gtf_gr[!(mcols(gtf_gr)$gene_id %in% filt_list)]

# save
rtracklayer::export.gff3(gtf_gr, str_replace(gtf_path, "\\.gff", ".filt.gff"))
  
# clean gene annotation
genes_info <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::select(gene_id, seqnames, start, end, strand) %>% 
  dplyr::left_join(., geneInfo %>% dplyr::select(gene_id, gene_name, gene_biotype, gene_description), by = "gene_id") %>% 
  dplyr::select(gene_id, seqnames, start, end, strand, gene_name, gene_biotype, gene_description)

# save as gene info
readr::write_csv(genes_info, str_replace(gtf_path, "\\.gff", ".filt.geneInfo.csv"))

# get exons
exons_gr <- 
  gtf_gr[mcols(gtf_gr)$type == "exon"] %>% 
  split(., mcols(.)$gene_id) %>% 
  GenomicRanges::reduce(.) %>% 
  unlist(.)
mcols(exons_gr)$gene_id <- names(exons_gr)
names(exons_gr) <- NULL

# save as reduced exons
exons_gr %>% 
  split(., mcols(.)$gene_id) %>% 
  saveRDS(., str_replace(gtf_path, "\\.gff", ".filt.reducedExons.RDS"))
