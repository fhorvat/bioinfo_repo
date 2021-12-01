### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/LTRs/potentially_young_LTRs")

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
library(stringdist)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# insertions table path
insertions_path <- file.path(inpath, "potentially_young_LTRs.edit_distance.10.csv")

# fasta files path
fasta_files_path <- list.files(inpath, ".*\\.fasta$")

######################################################## READ DATA
# read insertion table
insertions_tb <- readr::read_csv(insertions_path)

# read fasta files
insertions_seq <- 
  purrr::map(fasta_files_path, Biostrings::readDNAStringSet) %>% 
  do.call(c, .)

######################################################## MAIN CODE
# set rmsk as names of sequences
names(insertions_seq) <- str_remove(names(insertions_seq), ".*\\.")

# find sequences containing Ns
n_sequences <- 
  vmatchPattern(pattern = "N", subject = insertions_seq) %>% 
  unlist(.) %>% 
  names(.) %>% 
  unique(.)

# remove sequences containg N
insertions_seq <- insertions_seq[!names(insertions_seq) %in% n_sequences]

# find ORFs
# insertions_orfs <- predORF(insertions_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)
# saveRDS(insertions_orfs, file = file.path(outpath, str_c("potentially_young_LTRs.edit_distance.10", "ORFs.grl.RDS", sep = ".")))
insertions_orfs <- readRDS(file = file.path(outpath, str_c("potentially_young_LTRs.edit_distance.10", "ORFs.grl.RDS", sep = ".")))

# get tidy table
insertions_orfs_tb <- 
  insertions_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) 

# add info about ORFs to original table, filter
insertions_tb_full <- 
  insertions_tb %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::left_join(., insertions_orfs_tb, by = "rmsk_id") %>% 
  dplyr::arrange(desc(longest_orf_1 + longest_orf_2)) %>% 
  dplyr::select(insertion_coordinates:rmsk_id, longest_orf_1, longest_orf_2, everything(.))

# save 
readr::write_csv(insertions_tb_full, file.path(outpath, "potentially_young_LTRs.edit_distance.10.ORFs.csv"))

