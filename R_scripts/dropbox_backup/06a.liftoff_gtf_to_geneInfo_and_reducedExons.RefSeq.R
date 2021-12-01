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
gtf_path <- file.path(genome_dir, "annotation/Liftoff/MesAur1/RefSeq", "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

######################################################## MAIN CODE
# clean gene annotation
genes_info <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::mutate(gene_id = str_remove(ID, "gene-"), 
                gene_name = gene_id) %>%
  dplyr::select(gene_id, seqnames, start, end, strand, gene_name, gene_biotype, gene_description = description)

# save as gene info
readr::write_csv(genes_info, str_replace(gtf_path, "\\.gff", ".geneInfo.csv"))

# get exons
exons_gr <- 
  gtf_gr[mcols(gtf_gr)$type == "exon"] %>% 
  split(., mcols(.)$gene) %>% 
  GenomicRanges::reduce(.) %>% 
  unlist(.)
mcols(exons_gr)$gene_id <- names(exons_gr)
names(exons_gr) <- NULL

# save as reduced exons
exons_gr %>% 
  split(., mcols(.)$gene_id) %>% 
  saveRDS(., str_replace(gtf_path, "\\.gff", ".reducedExons.RDS"))