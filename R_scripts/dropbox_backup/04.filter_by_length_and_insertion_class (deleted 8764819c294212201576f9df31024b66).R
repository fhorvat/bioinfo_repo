### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/annotate_LINE1")

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

library(BSgenome.Mmusculus.UCSC.mm10)
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
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# LINE1s with ORF info path
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.annotated_exons.csv")

# L1Base path
l1base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

######################################################## READ DATA
# read LINE1 table
line1_tb <-
  readr::read_csv(line1_path) %>%
  dplyr::select_at(vars(-starts_with("X")))

# read L1Base
l1base_tb <- readr::read_csv(l1base_path)

######################################################## MAIN CODE
### add info about L1Base annotation overlap
# to GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.) 

# l1base to GRanges
l1base_gr <- GRanges(l1base_tb)

# find overlaps and percentage of overlap of each hit
line1_l1base_foverlaps <- findOverlaps(line1_gr, l1base_gr, ignore.strand = T)
# overlaps <- pintersect(line1_gr[queryHits(line1_l1base_foverlaps)], l1base_gr[subjectHits(line1_l1base_foverlaps)], ignore.strand = T)
# percentOverlap <- width(overlaps) / width(line1_gr[queryHits(line1_l1base_foverlaps)])

# extract ID's and genes
line1_l1base_tb <- 
  tibble(rmsk_id = line1_gr[queryHits(line1_l1base_foverlaps)]$rmsk_id,
         id = l1base_gr[subjectHits(line1_l1base_foverlaps)]$id) %>% 
  dplyr::group_by(rmsk_id) %>%
  dplyr::summarise(l1base_id_overlap = str_c(id, collapse = ", "))

# add info about overlaping exons to table
line1_tb %<>% dplyr::left_join(., line1_l1base_tb, by = "rmsk_id")

# save table with info about exon overlap
line1_tb %T>%
  readr::write_csv(., file.path(outpath, "LINE1.4000nt_plus.ORFs.annotated_exons_and_L1Base.csv"))

