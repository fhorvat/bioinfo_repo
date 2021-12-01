### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Piwil1_KO_analysis")

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

library(openxlsx)
library(rtracklayer)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

### genome files
# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# RefSeq gene info path
refseq_path <- file.path(genome_path, "annotation/Liftoff/MesAur1/RefSeq/hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.geneInfo.csv")

# ENSEMBL gene info path
ensembl_path <- file.path(genome_path, "annotation/Liftoff/MesAur1/ENSEMBL/ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.geneInfo.csv")

# Piwil3 results path
results_path.Piwil3 <- file.path(inpath, "Piwil3 DEGs.xlsx")

######################################################## READ DATA
# read RefSeq gene info
refseq_info <- readr::read_csv(refseq_path)

# read ENSEMBL gene info
ensembl_info <- readr::read_csv(ensembl_path)

# read results table 
results.Piwil3 <- openxlsx::read.xlsx(results_path.Piwil3) %>% as_tibble(.)

######################################################## MAIN CODE
# clean table
results.Piwil3_tidy <- 
  results.Piwil3 %>% 
  dplyr::select(gene_name = Gene.Name, lfc_piwil3_KO = M.Value)

# get GRanges
piwil3_gr <- 
  results.Piwil3_tidy %>% 
  dplyr::left_join(., refseq_info, by = "gene_name") %>% 
  GRanges(.)

# get ENSEMBL as GRanges
ensembl_gr <- 
  ensembl_info %>% 
  # dplyr::filter(gene_id %in% c("ENSMAUG00000012459",
  #                              "ENSMAUG00000013244",
  #                              "ENSMAUG00000014085",
  #                              "ENSMAUG00000019849",
  #                              "ENSMAUG00000007939",
  #                              "ENSMAUG00000001033")) %>%
  dplyr::mutate(gene_name = replace(gene_name, gene_name == "Shld2", "Fam35a"), 
                gene_name = replace(gene_name, gene_id == "ENSMAUG00000016142", "Usp24"),
                gene_name = replace(gene_name, gene_name == "Zfp551", "LOC101826595"), 
                gene_name = replace(gene_name, gene_id == "ENSMAUG00000013131", "Pou6f2")) %>% 
  GRanges(.) 

# overlap
overlap <- findOverlaps(piwil3_gr, ensembl_gr, ignore.strand = F)

# get as table
overlap_tb <- tibble(refseq_name = piwil3_gr[queryHits(overlap)]$gene_name,
                     refseq_seqnames = as.character(seqnames(piwil3_gr[queryHits(overlap)])), 
                     refseq_start = start(piwil3_gr[queryHits(overlap)]), 
                     refseq_end = end(piwil3_gr[queryHits(overlap)]), 
                     
                     ensembl_seqnames = as.character(seqnames(ensembl_gr[subjectHits(overlap)])), 
                     ensembl_start = start(ensembl_gr[subjectHits(overlap)]), 
                     ensembl_end = end(ensembl_gr[subjectHits(overlap)]), 
                     
                     gene_id = ensembl_gr[subjectHits(overlap)]$gene_id, 
                     gene_name = ensembl_gr[subjectHits(overlap)]$gene_name) %>% 
  dplyr::filter(refseq_name != gene_name)






