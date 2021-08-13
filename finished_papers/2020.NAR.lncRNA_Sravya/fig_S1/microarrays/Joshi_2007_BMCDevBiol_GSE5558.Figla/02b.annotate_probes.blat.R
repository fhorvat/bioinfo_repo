### INFO: BLATs probe sequences to the genome, gets overlaps with annotated genes
### DATE: Sat Aug 17 18:47:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Joshi_2007_BMCDevBiol_GSE5558/align_probes/blat")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(GEOquery)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# blat .psl path
blat_psl_path <- list.files(path = inpath, pattern = ".*\\.score\\.psl$", full.names = T)

# blat .bed path (with strand info)
blat_bam_path <- list.files(path = inpath, pattern = ".*\\.bam$", full.names = T)

# tidy annotation table path
annotation_tidy_path <- file.path(inpath, "../../GPL3771.annot.tidy.csv")

# reduced exon coordinates path
exons_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS"

# gene info path
gene_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv"

######################################################## READ DATA
# read blat .psl
blat_psl <- 
  readr::read_delim(file = blat_psl_path, delim = "\t", col_names = c("seqnames", "start", "end", "subject_coords", "score", "perc_identity")) %>% 
  dplyr::mutate(probe = str_remove(subject_coords, ":.*"), 
                start = start + 1) %>% 
  dplyr::select(seqnames, start, end, probe, score)

# read blat .bed
blat_bam <- GenomicAlignments::readGAlignments(blat_bam_path, use.names = T)

# read annotation table
annotation_tidy <- readr::read_csv(annotation_tidy_path)

# read exon coordinates
exons_gr <- readRDS(exons_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
### tidy probe alignments
# check if .psl and .bam files are parallel
if(all(names(blat_bam) == blat_psl$probe) & all(seqnames(blat_bam) == blat_psl$seqnames) & 
   all(start(blat_bam) == blat_psl$start) & all(end(blat_bam) == blat_psl$end)){
  
  # add score to alignment file
  mcols(blat_bam)$score <- blat_psl$score
  
}else{
  
  # stop with error message
  stop("Psl and bam file from BLAT are not parallel!")
  
}

# set unique names to alignments
names(blat_bam) <- make.unique(names(blat_bam))

# for each probe get top score
probe_top_scores <- 
  tibble(probe = names(blat_bam), 
         score = mcols(blat_bam)$score) %>% 
  dplyr::mutate(probe_original = str_remove(probe, "\\.[0-9]+")) %>% 
  dplyr::group_by(probe_original) %>% 
  dplyr::filter(score == max(score)) %>% 
  dplyr::ungroup(.)

# filter top scoring alignments for each probe, get ranges of only aligned part
probe_gr <- 
  blat_bam[names(blat_bam) %in% probe_top_scores$probe] %>% 
  GenomicRanges::grglist(.)

# find overlaps between probes and alignments
probe_gene_overlaps <- findOverlaps(probe_gr, exons_gr, ignore.strand = T)

# extract probes and gene info
annotation_blat <- 
  tibble(GenBank = probe_gr[queryHits(probe_gene_overlaps)] %>% names,
         gene_id = exons_gr[subjectHits(probe_gene_overlaps)] %>% names(.)) %>% 
  dplyr::left_join(., gene_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%  
  dplyr::mutate(GenBank = str_remove(GenBank, "\\.[0-9]+$")) %>% 
  dplyr::group_by(GenBank) %>% 
  dplyr::summarise(gene_name = str_c(gene_name, collapse = "."), 
                   gene_id = str_c(gene_id, collapse = ".")) %>% 
  dplyr::ungroup(.)

# add annotation by BLAT to table
annotation_tb <- 
  annotation_tidy %>% 
  dplyr::rename(gene_description = gene_name) %>% 
  dplyr::left_join(., annotation_blat, by = "GenBank")


