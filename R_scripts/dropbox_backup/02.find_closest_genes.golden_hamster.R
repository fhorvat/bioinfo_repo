### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/MYSERV/FLI/golden_hamster.Siomi")

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
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
 
######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_name <- basename(inpath)
genome_dir <- file.path(inpath, str_c("../../genomes/", genome_name))

# # joined repeatMasker path
# rmsk_path <- list.files(genome_dir, "rmsk\\..*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)
# 
# # clean repeatMasker path 
# rmsk_clean_path <- list.files(genome_dir, "rmsk\\..*\\.clean\\.fa\\.out\\.gz", full.names = T)

# top 10 MYSERV coordinates
myserv_bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/MYSERV6/annotation"
myserv_bed_path <- file.path(myserv_bed_path, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.bed")

# gene annotation .gtf path
gtf_path <- list.files(genome_dir, "ensembl\\.99.*\\.UCSCseqnames\\.gff", full.names = T)

# gene info path
genes_info_path <- list.files(genome_dir, pattern = "ensembl\\.99.*\\.UCSCseqnames.geneInfo.csv$", full.names = T)

# chinese hamster gene info
chamster_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"
chamster_path <- list.files(chamster_path, pattern = "ensembl\\.99.*\\.UCSCseqnames.geneInfo.csv$", full.names = T)
  
######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")
# 
# # read clean repeatMasker
# rmsk_clean <- readr::read_delim(rmsk_clean_path, delim  = "\t")

# read MYSERV coordinates
myserv_gr <- rtracklayer::import(myserv_bed_path)

# read .gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read chinese hamster genes info
chamster_genes_info <- readr::read_csv(chamster_path)

######################################################## MAIN CODE
### clean data
# get only genes
gene_gr <- gtf_gr[mcols(gtf_gr)$type == "gene"]
gene_gr <- gene_gr[mcols(gene_gr)$gene_biotype == "protein_coding"]

# for each insertion find closest up-stream and downstream gene
mcols(myserv_gr)$preceded <- mcols(gene_gr)$gene_id[GenomicRanges::precede(myserv_gr, gene_gr, ignore.strand = T)]
mcols(myserv_gr)$followed <- mcols(gene_gr)$gene_id[GenomicRanges::follow(myserv_gr, gene_gr, ignore.strand = T)]

# add info about genes
myserv_gene_tb <- 
  myserv_gr %>% 
  as_tibble(.) %>% 
  tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ") %>% 
  dplyr::select(coordinates, name, followed, preceded) %>% 
  tidyr::pivot_longer(cols = -c(coordinates, name), names_to = "position", values_to = "gene_id") %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% 
  dplyr::filter(!is.na(gene_name)) %>% 
  dplyr::left_join(., chamster_genes_info %>% 
                     tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ") %>% 
                     dplyr::select(gene_id_chamster = gene_id, gene_name, coordinates_chamster = coordinates),
                   by = "gene_name") %>% 
  tidyr::pivot_wider(id_cols = c(name, coordinates), names_from = position, values_from = gene_id_chamster)


