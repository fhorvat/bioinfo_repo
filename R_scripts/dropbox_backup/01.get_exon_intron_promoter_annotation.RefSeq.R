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

# gtf path
gtf_path <- file.path(genome_dir, "annotation/Liftoff/MesAur1/RefSeq", "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff")

# repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, "\t")

######################################################## MAIN CODE
# clean gene annotation
genes_info <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::mutate(gene_id = str_remove(ID, "gene-")) %>% 
  dplyr::select(gene_id, seqnames, start, end, strand, gene_name = gene_id, gene_biotype, gene_description = description)

# get exons
exons_gr <- 
  gtf_gr[mcols(gtf_gr)$type == "exon"] %>% 
  split(., mcols(.)$gene) %>% 
  GenomicRanges::reduce(.) %>% 
  unlist(.)
mcols(exons_gr)$gene_id <- names(exons_gr)
names(exons_gr) <- NULL

# get introns
introns_gr <- endoapply(exons_gr, function(x) x %>% gaps(.) %>% .[start(.) != 1])
introns_gr <- unlist(introns_gr)
mcols(introns_gr)$gene_id <- names(introns_gr)
names(introns_gr) <- NULL

# get promoters (0.5 kb)
promoters_500bp_gr <- 
  genes_info %>% 
  GRanges(.) %>% 
  GenomicRanges::promoters(., upstream = 500, downstream = 0)
mcols(promoters_500bp_gr) <- mcols(promoters_500bp_gr)$gene_name
names(mcols(promoters_500bp_gr)) <- "gene_id"

# get promoters (1 kb)
promoters_1000bp_gr <- 
  genes_info %>% 
  GRanges(.) %>% 
  GenomicRanges::promoters(., upstream = 1000, downstream = 0)
mcols(promoters_1000bp_gr) <- mcols(promoters_1000bp_gr)$gene_name
names(mcols(promoters_1000bp_gr)) <- "gene_id"

# save as list of GRanges
list(exons = exons_gr, introns = introns_gr, promoters_500bp = promoters_500bp_gr, promoters_1000bp = promoters_1000bp_gr) %>% 
  saveRDS(., file = file.path(outpath, "RefSeq_MesAur1.exon_intron_promoters.coords.GRanges.RDS"))
