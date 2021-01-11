### INFO: 
### DATE: Mon Jul 08 21:52:06 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/MesAur1.LINE1_ORF1/expand_150b_around_ORF1")

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

library(BSgenome.Maur.UCSC.MesAur1)
library(Biostrings)
library(seqinr)
library(msa)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# blast results path
blast_results_path <- list.files(path = inpath, pattern = ".*\\.blastn\\.txt", full.names = T)

######################################################## READ DATA
# read blast results
blast_results <- purrr::map(blast_results_path, function(path){
  
  # read table and name it
  readr::read_delim(file = path, delim = "\t", col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                             "mismatches", "gap_open", 
                                                             "query_start", "query_end", "subject_start", "subject_end", 
                                                             "e_value", "bit_score")) %>% 
    arrange(-alignment_length) 
  
}) %>% 
  magrittr::set_names(., basename(blast_results_path) %>% str_remove(., "\\.blastn\\.txt"))

######################################################## MAIN CODE
# get sequences, MSA and phylotree
blast_seq <- 
  blast_results[[1]] %>% 
  dplyr::filter(alignment_length > 250) %>%
  dplyr::arrange(desc(identity_perc)) %>% 
  dplyr::select(seqnames = subject_id, subject_start, subject_end, query = query_id) %>% 
  dplyr::mutate(strand = ifelse(subject_start < subject_end, "+", "-")) %>% 
  dplyr::mutate(start = ifelse(strand == "+", subject_start, subject_end), 
                end = ifelse(strand == "+", subject_end, subject_start)) %>% 
  dplyr::select(seqnames, start, end, strand, query)

# to GRanges
line1_gr <-
  blast_seq %>%
  GRanges(.)
names(line1_gr) <- str_c(seqnames(line1_gr), ":", start(line1_gr), "-", end(line1_gr))
seqlevels(line1_gr) <- seqlevels(BSgenome.Maur.UCSC.MesAur1)
seqlengths(line1_gr) <- seqlengths(BSgenome.Maur.UCSC.MesAur1)
line1_gr <- trim(line1_gr)
names(line1_gr) <- str_replace_all(names(line1_gr), ":|-", " ")

# get sequences
line1_seq <- getSeq(BSgenome.Maur.UCSC.MesAur1, line1_gr)

# save sequence
Biostrings::writeXStringSet(line1_seq, file = file.path(outpath, str_c("LINE1.CriGriPICR.ORF1_expand_150b.L1-1_CGr.consensus", "BLAST_hits_in_MesAur1", "fasta", sep = ".")))


### MSA together
# msa
seq_msa <- msa::msa(line1_seq, method = "ClustalOmega")

# # write as fasta
# seq_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("LINE1.CriGriPICR.ORF1_expand_150b.L1-1_CGr.consensus", "BLAST_hits_in_MesAur1", "ClustalO.msa.fasta", sep = ".")))

# convert to seqinr format
msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")

# calculate distance matrix
msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")

# get tree and plot
msa_tree <- ape::nj(msa_dist)
msa_tree <- ape::ladderize(msa_tree)

# save with labels
pdf(file = file.path(outpath, 
                     str_c("LINE1.CriGriPICR.ORF1_expand_150b.L1-1_CGr.consensus", "BLAST_hits_in_MesAur1", "ClustalO.msa", "phylo_tree", "pdf", sep = ".")), 
    width = 20, height = 60)
plot.phylo(msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
# edgelabels(round(msa_tree$edge.length, 3), bg = "black", col = "white", font = 2)
axisPhylo(backward = F)
dev.off()

# table
line1_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(-width) %>% 
  tidyr::unite(col = coordinates, "seqnames", "start", "end", sep = " ") %>% 
  dplyr::right_join(., tibble(coordinates = msa_tree$tip.label), by = "coordinates") %T>% 
  readr::write_csv(., path = file.path(outpath, str_c("LINE1.CriGriPICR.ORF1_expand_150b.L1-1_CGr.consensus", "BLAST_hits_in_MesAur1", "csv", sep = ".")))



