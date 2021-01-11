### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq/Data/Mapped/STAR_Siomi/7_filter_reads_by_coordinates")

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
library(GenomicRanges)
library(GenomicAlignments)

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

# bam path 
bam_path <- file.path(inpath, "s_testis_Mov10l_WT_8.5_dpp.24to31nt.merged.bam")

######################################################## READ DATA
# read gtf
gtf_gr <- rtracklayer::import(gtf_path)

# read bam
bam_gr <- readGAlignmentsList(bam_path, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = NA)))

######################################################## MAIN CODE
### get coding regions
# get exons
exons_gr <- 
  gtf_gr[mcols(gtf_gr)$type == "exon"] %>% 
  split(., mcols(.)$gene_id)
# exons_gr <- exons_gr[names(exons_gr) == "ENSMAUG00000013712"]
exons_gr <- unlist(exons_gr)

### get read coordinates
# unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
bam_gr <- unlist(bam_gr)
names(bam_gr) <- mcols(bam_gr)$qname
bam_gr <-
  GenomicRanges::grglist(bam_gr) %>%
  unlist(.)

### find reads origin
# find overlaps
overlaps <- findOverlaps(bam_gr, exons_gr, ignore.strand = T)

# merge reads to overlap
clusters_gr <- GenomicRanges::reduce(bam_gr, min.gapwidth = 100, ignore.strand = F)

# count reads per cluster
read_counts <- countOverlaps(clusters_gr, bam_gr, ignore.strand = F)

# add counts to clusters
clusters_tb <- 
  clusters_gr %>% 
  as_tibble(.) %>% 
  dplyr::mutate(counts = read_counts) %>% 
  dplyr::arrange(-counts)



