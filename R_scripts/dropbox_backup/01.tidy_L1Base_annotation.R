### INFO: 
### DATE: Wed Oct 09 15:44:16 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1Base")

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


## DOCUMENTATION
# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1Base"

# LINE1s coordinates
line1_path <- file.path(documentation_path, "2019_RNAi_piRNA_file3_L1_table.csv")


######################################################## READ DATA
# read LINE1 table
line1_tb <-
  readr::read_csv(line1_path) %>%
  dplyr::select_at(vars(-starts_with("X")))

######################################################## MAIN CODE
# create GRanges
line1_gr <-
  line1_tb %>%
  dplyr::rename(rmsk_id = L1base_id,
                seqnames = chrName, start = chrStart, end = chrEnd) %>%
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>%
  GRanges(.)
names(line1_gr) <- mcols(line1_gr)$rmsk_id

# save as .RDS
# saveRDS(line1_gr, file.path(documentation_path, "LINE1_annotation.L1Base.GRanges.RDS"))

# # save as .bed
# rtracklayer::export.bed(line1_gr, file.path(documentation_path, "LINE1_annotation.L1Base.bed"))

# save as SAF
line1_saf <- 
  line1_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
  readr::write_delim(., file.path(outpath, "LINE1_annotation.L1Base.saf"), delim = "\t")
