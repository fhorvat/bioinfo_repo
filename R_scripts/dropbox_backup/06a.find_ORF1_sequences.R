### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/chinese_hamster/CriGri_PICR.GCA_003668045/LINE1")

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

library(BSgenome.Cgriseus.Ensembl.CriGriPICR)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LINE1s ORFs path
line1_orfs_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.grl.RDS")

# LINE1 table path
line1_tb_path <- file.path(inpath, "LINE1.5K_to_7k.ORFs.annotated_exons.csv")

######################################################## READ DATA
# read LINE1 ORFs
line1_orfs <- readRDS(line1_orfs_path)

# read LINE1 table
line1_tb <- 
  readr::read_csv(line1_tb_path) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

######################################################## MAIN CODE
# get coordinates of the 2 longest ORFs for each LINE1 
# find ORF1 coordinates
# expand 1.5 kb upstream
line1_orfs_tb <- 
  line1_orfs[names(line1_orfs) %in% line1_tb$rmsk_id] %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, orf_start = start, orf_end = end, orf_width = width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>%
  dplyr::top_n(2, orf_width) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(orf_width < 2000) %>% 
  dplyr::left_join(., line1_tb %>% dplyr::select(seqnames, start, end, strand, rmsk_id, repName), by = "rmsk_id") %>% 
  dplyr::select(seqnames, start, end, strand, orf_start, orf_end, orf_width, repName, rmsk_id) %>% 
  dplyr::mutate(orf_abs_start = ifelse(strand == "+", start + orf_start - 1, end - orf_start + 1), 
                orf_abs_end = ifelse(strand == "+", orf_abs_start + orf_width - 1, orf_abs_start - orf_width + 1), 
                seq_start = ifelse(strand == "+", orf_abs_start, orf_abs_end), 
                seq_end = ifelse(strand == "+", orf_abs_end, orf_abs_start)) %>% 
  dplyr::select(seqnames, start = seq_start, end = seq_end, strand, orf_width, rmsk_id, repName)

# to GRanges
line1_orfs_gr <-
  line1_orfs_tb %>%
  GRanges(.)
names(line1_orfs_gr) <- str_c(mcols(line1_orfs_gr)$repName, mcols(line1_orfs_gr)$rmsk_id, sep = ".")
names(line1_orfs_gr) <- str_replace_all(names(line1_orfs_gr), "-", "_")
seqlevels(line1_orfs_gr) <- seqlevels(BSgenome.Cgriseus.Ensembl.CriGriPICR)
seqlengths(line1_orfs_gr) <- seqlengths(BSgenome.Cgriseus.Ensembl.CriGriPICR)
line1_orfs_gr <- trim(line1_orfs_gr)
mcols(line1_orfs_gr)$repName <- str_replace(mcols(line1_orfs_gr)$repName, "L1-2", "L1-1")

# get sequences
line1_orfs_seq <- getSeq(BSgenome.Cgriseus.Ensembl.CriGriPICR, line1_orfs_gr)
# line1_orfs_seq <- line1_orfs_seq[order(names(line1_orfs_seq))]

# save as fasta
Biostrings::writeXStringSet(x = line1_orfs_seq, filepath = file.path(outpath, "LINE1.CriGriPICR.ORF1.fasta"))


### get consensus by repeat names
# split by repeat names
line1_seq_list <- split(line1_orfs_seq, line1_orfs_gr$repName)

## MSA and get consesus
line1_seq_consensus <- purrr::map(names(line1_seq_list), function(seq_name){
  
  # do the MSA
  seq_msa <- msa::msa(line1_seq_list[[seq_name]], method = "ClustalOmega")
  
  # get the consensus
  seq_consensus <- 
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>% 
    toupper(.) %>% 
    Biostrings::DNAStringSet(.)
  names(seq_consensus) <- str_c(seq_name, ".consensus")
  
  # write as fasta
  Biostrings::writeXStringSet(seq_consensus, file = file.path(outpath, str_c("LINE1.CriGriPICR.ORF1", seq_name, "consensus.fasta", sep = ".")))
  
  # return
  return(seq_consensus)
  
})
