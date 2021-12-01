### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/Papd7_KO/polyA_tails")

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
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# original annotation table
utr_path <- file.path(inpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.3pUTR.annotation.csv")

# bam path
bam_path <- file.path(inpath, "bam_subset")

# bam list
bam_list <- list.files(bam_path, pattern = ".*\\.bam$", full.names = T)

######################################################## READ DATA
# read UTR coordinates
utr_tb <- readr::read_csv(utr_path)

# read bam
bam_gr <- GenomicAlignments::readGAlignmentsList(bam_list[1], 
                                                 param = ScanBamParam(what = c("qname", "seq")))

######################################################## MAIN CODE
# get UTR annotations as GRanges
utr_gr <- GRanges(utr_tb)

# get soft-clipped reads
bam_filt <- 
  bam_gr %>% 
  unlist(.) %>% 
  .[str_detect(cigar(.), "35S")]

# overlap with UTR annotation
overlap <- findOverlaps(bam_filt, utr_gr, ignore.strand = T)

# table
bam_tb <- 
  tibble(seqnames = as.character(seqnames(bam_filt[queryHits(overlap)])), 
         start = as.numeric(start(bam_filt[queryHits(overlap)])), 
         end = as.numeric(end(bam_filt[queryHits(overlap)]))) %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::mutate(seqnames = as.character(seqnames(utr_gr[subjectHits(overlap)])), 
                start = as.numeric(start(utr_gr[subjectHits(overlap)])), 
                end = as.numeric(end(utr_gr[subjectHits(overlap)]))) %>% 
  tidyr::unite(col = coordinates_utr, seqnames, start, end, sep = " ")
  
  
  