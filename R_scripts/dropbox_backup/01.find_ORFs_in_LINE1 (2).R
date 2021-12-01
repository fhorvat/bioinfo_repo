### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/chinese_hamster/LINE1")

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

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/chinese_hamster/criGriChoV2.CHOK1GS_HDv1.GCA_900186095.1"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.criGriChoV2.20200423.joined_rmsk_id.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# get widths
rmsk_tb %<>%    
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "+", strand)) %>% 
  GRanges(.) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

# filter LINE1s
line1_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repFamily == "L1", 
                width > 4000)

# to GRanges
line1_gr <- GRanges(line1_tb) 

# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Cgriseus.UCSC.criGriChoV2, names = line1_gr)
names(line1_seq) <- line1_gr$rmsk_id

# find sequences containing Ns
n_sequences <- 
  vmatchPattern(pattern = "N", subject = line1_seq) %>% 
  unlist(.) %>% 
  names(.) %>% 
  unique(.)

# remove sequences containg N
line1_seq <- line1_seq[!names(line1_seq) %in% n_sequences]

# find ORFs
line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)
saveRDS(line1_orfs, file = file.path(outpath, "LINE1.4000nt_plus.ORFs.grl.RDS"))

# get tidy table
line1_orfs_tb<- 
  line1_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) 

# add info about ORFs to original table, filter
line1_tb %>% 
  dplyr::left_join(., line1_orfs_tb, by = "rmsk_id") %>% 
  dplyr::filter(rmsk_id %in% names(line1_seq)) %T>%
  readr::write_csv(., file.path(outpath, "LINE1.4000nt_plus.ORFs.csv"))

