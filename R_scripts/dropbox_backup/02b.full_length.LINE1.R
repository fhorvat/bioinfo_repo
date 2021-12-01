### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1_joined_rmskID/mesAur1.L1_joined_rmskID_20190929.FH")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.MesAur1.0.20160612.RefSeq.joined_rmsk_ID.fa.out.gz")

# reduced exons path
exons_path <- file.path(genome_dir, "ensembl.93.MesAur1.0.20180920.reducedExons.RDS")

# fasta index path
fai_path <- file.path(genome_dir, "MesAur1.0.fa.gz.fai")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# get exon ranges as GRanges
exons_gr <- unlist(exons_gr)
# seqlevels(exons_gr) <- fai_tb$seqnames

# repeatMasker to GRanges, get only LINE1's
line1_whole <- 
  rmsk_tb %>%    
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.) %>% 
  as_tibble(.) %>% 
  dplyr::filter(repFamily == "L1", 
                insertion_class == "whole", 
                width > 6000) %>% 
  dplyr::rename(rmsk_id = rmsk_ID)

# export bed
line1_gr <- 
  line1_whole %>%
  GRanges(.) %T>%
  rtracklayer::export.bed(., file.path(outpath, "LINE1_whole.joined_rmsk_ID.bed"))


### find whole LINE1s which overlap exons of annotated genes
# overlap
line1_exon_foverlaps <- GenomicRanges::findOverlaps(line1_gr, exons_gr, ignore.strand = T, minoverlap = 1)

# extrack ID's and genes
line1_exon_tb <- 
  tibble(rmsk_id = line1_gr[queryHits(line1_exon_foverlaps)]$rmsk_id, 
         overlaps_ensembl_id = exons_gr[subjectHits(line1_exon_foverlaps)]$gene_id)

# export csv
line1_whole %>%
  dplyr::left_join(., line1_exon_tb, by = "rmsk_id") %T>% 
  readr::write_csv(., file.path(outpath, "LINE1_whole.joined_rmsk_ID.csv"))


