### INFO: 
### DATE: Fri Nov 20 10:05:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/small_RNAseq")

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
# set outpath
outpath <- getwd()

# set inpath
inpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# ENSEMBL path
features_ensembl <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")

# cluster table path
cluster_tb_list <- list.files(inpath, "piRNA_clusters.*\\.20210513\\.csv", full.names = T)

######################################################## READ DATA
# read .gtf
gtf_gr <- rtracklayer::import(features_ensembl)

######################################################## MAIN CODE
# get miRNA annotation
mirna_gr <- gtf_gr[mcols(gtf_gr)$gene_biotype == "miRNA" & mcols(gtf_gr)$type == "gene"]

### loop through tables
purrr::map(cluster_tb_list, function(cluster_tb_path){
  
  # read cluster table
  cluster_tb <- readr::read_csv(cluster_tb_path)

  # get cluster coordinates
  cluster_gr <- 
    cluster_tb %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
    dplyr::mutate(strand = "*") %>% 
    GRanges(.)
  
  # overlap clusters with miRNA
  overlaps <- findOverlaps(cluster_gr, mirna_gr, ignore.strand = T)
  
  # get coordinates of overlaping clusters
  cluster_coordinates <- 
    mcols(cluster_gr[queryHits(overlaps)])$coordinates %>% 
    unique(.) 
  
  # add annotation to table
  cluster_tb %<>% 
    dplyr::mutate(miRNA = ifelse(coordinates %in% cluster_coordinates, "yes", "no")) %T>% 
    readr::write_csv(., file = str_replace(cluster_tb_path, "\\.csv$", ".annotated_miRNA.csv"))
  
})
