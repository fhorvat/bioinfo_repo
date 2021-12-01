### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/rmsk")

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

# LINE1s rmsk path
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.annotated_exons_and_L1Base.csv")

######################################################## READ DATA
# read LINE1 table
line1_tb <- readr::read_csv(line1_path)

######################################################## MAIN CODE
# filter table and save
line1_out <- 
  line1_tb %>% 
  dplyr::filter(width >= 5000, width <= 7000, 
                insertion_class %in% c("whole", "within")) %>% 
  dplyr::select(-c(repClass, repFamily)) %T>%
  readr::write_csv(., file.path(outpath, "LINE1.LINE1.5K_to_7k.ORFs.annotated_exons_and_L1Base.csv"))

# to GRanges
line1_gr <- 
  line1_out %>% 
  GRanges(.)
names(line1_gr) <- mcols(line1_gr)$rmsk_id

# # save
# saveRDS(line1_gr, file = file.path(outpath, str_c("LINE1_annotation.rmsk.GRanges.RDS")))

# # save as .bed
# rtracklayer::export.bed(line1_gr, file.path(outpath, "LINE1_annotation.rmsk.bed"))

# save as SAF
line1_saf <- 
  line1_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
  readr::write_delim(., file.path(outpath, "LINE1_annotation.rmsk.saf"), delim = "\t")

