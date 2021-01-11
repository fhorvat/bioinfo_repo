### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/methylation/non_repetitive")

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
genome_dir <- file.path(genome_dir, "annotation/Liftoff/MesAur1/ENSEMBL")

# gtf path
gtf_path <- file.path(genome_dir, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.gff")

# geneInfo path
geneInfo_path <- file.path(genome_dir, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.geneInfo.csv")

# get gene FPKM path
fpkm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression.Siomi"
fpkm_path <- file.path(fpkm_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.FPKM_mean.csv")

# repeatMasker path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read genesInfo
genes_info <- readr::read_csv(geneInfo_path)
  
# read FPKM
fpkm_tb <- readr::read_csv(fpkm_path)

# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
### clean data
# get list of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter repeatMasker, get GRanges
rmsk_gr <- 
  rmsk_tb %>% 
  # dplyr::filter(repClass %in% c("LINE", "LTR", "SINE", "DNA", "Retroposon")) %>%
  dplyr::filter(!(repClass %in% c("Low_complexity", "Simple_repeat"))) %>%
  dplyr::mutate(repClass = replace(repClass, !(repClass %in% c("LINE", "LTR", "SINE")), "other_repeat")) %>% 
  dplyr::mutate(gene_id = rmsk_id) %>% 
  dplyr::select(seqnames, start, end, strand, gene_id, repClass) %>% 
  GRanges(.)

# annotate genes based on expression
genes_fpkm <- 
  fpkm_tb %>% 
  dplyr::mutate(expression = ifelse(Mov10l_WT >= 0.5, "expressed", "non_expressed")) %>% 
  dplyr::select(gene_id, expression)


### get coding regions
# get exons
exons_gr <- 
  gtf_gr[mcols(gtf_gr)$type == "exon"] %>% 
  split(., mcols(.)$gene_id) %>% 
  GenomicRanges::reduce(.) 
exons_gr <- exons_gr[names(exons_gr) %in% protein_genes]

# get introns
introns_gr <- endoapply(exons_gr, function(x) x %>% gaps(.) %>% .[start(.) != 1])
introns_gr <- unlist(introns_gr)
mcols(introns_gr)$gene_id <- names(introns_gr)
names(introns_gr) <- NULL

# unlist exons
exons_gr <- unlist(exons_gr)
mcols(exons_gr)$gene_id <- names(exons_gr)
names(exons_gr) <- NULL


### get promoters
# get promoters (1 kb)
promoters_gr <- 
  genes_info %>% 
  dplyr::filter(gene_id %in% protein_genes) %>% 
  dplyr::left_join(., genes_fpkm, by = "gene_id") %>% 
  dplyr::filter(!is.na(expression)) %>% 
  GRanges(.) %>% 
  GenomicRanges::promoters(., upstream = 1000, downstream = 0)
mcols(promoters_gr) <- mcols(promoters_gr)[, c("gene_id", "expression")]

# save as list of GRanges
list(LINE = rmsk_gr[rmsk_gr$repClass == "LINE"], 
     LTR = rmsk_gr[rmsk_gr$repClass == "LTR"],
     SINE = rmsk_gr[rmsk_gr$repClass == "SINE"], 
     rmsk_other = rmsk_gr[rmsk_gr$repClass == "other_repeat"],
     promoters_expressed = promoters_gr[mcols(promoters_gr)$expression == "expressed"], 
     promoters_non_expressed = promoters_gr[mcols(promoters_gr)$expression == "non_expressed"], 
     exons = exons_gr, 
     introns = introns_gr) %>% 
  GRangesList(.) %>% 
  saveRDS(., file = file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.exons_introns_promoters.RDS"))
