### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/chinese_hamster/criGriChoV2.CHOK1GS_HDv1.GCA_900186095.1/LINE1")

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

library(BSgenome.Cgriseus.UCSC.criGriChoV2)
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

# genome path
genome_dir <- "/common/DB/genome_reference/chinese_hamster/criGriChoV2.CHOK1GS_HDv1.GCA_900186095.1"

# repeatMasker
rmsk_path <- file.path(genome_dir, "rmsk.criGriChoV2.20200423.clean.fa.out.gz")

# LINE1s rmsk path
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.annotated_exons.csv")

######################################################## READ DATA
# read LINE1 table
line1_tb <- readr::read_csv(line1_path)

# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# filter table and save
line1_out <- 
  line1_tb %>% 
  dplyr::filter(width >= 5000, width <= 7000, 
                insertion_class %in% c("whole", "within")) %>% 
  dplyr::arrange(-longest_orf_1) %>% 
  dplyr::select(-c(repClass, repFamily)) %T>%
  readr::write_csv(., file.path(outpath, "LINE1.5K_to_7k.ORFs.annotated_exons.csv"))

# ORF1 at least 370 amino-acids long and ORF2 should be at least 1200 amino acids long

# # to GRanges
# line1_gr <- 
#   line1_out %>% 
#   GRanges(.)
# names(line1_gr) <- mcols(line1_gr)$rmsk_id
# 
# # # save
# # saveRDS(line1_gr, file = file.path(outpath, str_c("LINE1_annotation.rmsk.GRanges.RDS")))
# 
# # # save as .bed
# # rtracklayer::export.bed(line1_gr, file.path(outpath, "LINE1_annotation.rmsk.bed"))
# 
# # save as SAF
# line1_saf <- 
#   line1_gr %>% 
#   as_tibble(.) %>% 
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
#   readr::write_delim(., file.path(outpath, "LINE1_annotation.rmsk.saf"), delim = "\t")
# 
#  