### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set parameters and paths 
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression.joined_rmsk_ID"
genome_version <- "mm10.L1_joined_rmskID_20190929.FH"
experiment <- "Dicer_Mili_KO"
single_end <- FALSE

# set working dir
setwd(file.path(base_path, "Documentation", genome_version))

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
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz")

# reduced exons path
exons_path <- file.path(genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# repeatMasker to GRanges, get only LINE1's
line1_whole <- 
  rmsk_tb %>%    
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.) %>% 
  as_tibble(.) %>% 
  dplyr::filter(repFamily == "L1", 
                insertion_class == "whole", 
                width > 6000)

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
         overlaps_ensembl_id = names(exons_gr[subjectHits(line1_exon_foverlaps)])) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(overlaps_ensembl_id = str_c(overlaps_ensembl_id, collapse = ", "))

# export csv
line1_whole_tb <-
  line1_whole %>%
  dplyr::left_join(., line1_exon_tb, by = "rmsk_id") %T>% 
  readr::write_csv(., file.path(outpath, "LINE1_whole.joined_rmsk_ID.csv"))


