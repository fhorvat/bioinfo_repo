### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/chinese_hamster/CriGri_PICR.GCA_003668045/LINE1")

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
library(DESeq2)

library(BSgenome.Cgriseus.Ensembl.CriGriPICR)
library(seqinr)
library(Biostrings)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/chinese_hamster/CriGri_PICR.GCA_003668045"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.CriGri_PICR.20200428.joined_rmsk_id.fa.out.gz")

# reduced exons path
exons_path <- file.path(genome_dir, "ensembl.99.CriGri_PICR.20200428.UCSCseqnames.reducedExons.RDS")

# gene info path
genes_info_path <- file.path(genome_dir, "ensembl.99.CriGri_PICR.20200428.UCSCseqnames.geneInfo.csv")

# LINE1s with ORF info path
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.csv")

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read LINE1 table
line1_tb <- 
  readr::read_csv(line1_path) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

######################################################## MAIN CODE
### add info about Ensembl exon overlap
# to GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.) 

# find overlaps with exons of annotated genes
line1_exon_foverlaps <- GenomicRanges::findOverlaps(line1_gr, exons_gr, ignore.strand = T)

# extract ID's and genes
line1_exon_tb <- 
  tibble(rmsk_id = as.character(line1_gr[queryHits(line1_exon_foverlaps)]$rmsk_id), 
         ensembl_id_overlap = names(exons_gr[subjectHits(line1_exon_foverlaps)])) %>% 
  dplyr::left_join(., genes_info, by = c("ensembl_id_overlap" = "gene_id")) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(ensembl_id_overlap = str_c(ensembl_id_overlap, collapse = ", "), 
                   ensembl_gene_biotype_overlap = str_c(gene_biotype, collapse = ", "))

# add info about overlaping exons to table
line1_tb %<>% dplyr::left_join(., line1_exon_tb, by = "rmsk_id")

# save table with info about exon overlap
line1_tb %T>%
  readr::write_csv(., file.path(outpath, "LINE1.4000nt_plus.ORFs.annotated_exons.csv"))

  

