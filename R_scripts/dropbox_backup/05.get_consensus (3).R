### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV.repeatMasker/annotation")

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

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)
library(msa)
library(DECIPHER)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# insertions table path
insertions_path <- file.path(inpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs 201213.xlsx")

######################################################## READ DATA
# read insertion table
insertions_tb <- openxlsx::read.xlsx(insertions_path) %>% as_tibble(.)

######################################################## MAIN CODE
# tidy table
insertions_tb %<>% 
  dplyr::slice(1:10) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(start = as.numeric(start), 
                end = as.numeric(end), 
                width = end - start + 1)

# to GRanges
insertions_gr <-
  insertions_tb %>%
  GRanges(.)
names(insertions_gr) <- str_c("MYSERV", mcols(insertions_gr)$rmsk_id, sep = ".")

# save as .bed
rtracklayer::export.bed(insertions_gr, file.path(outpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.bed"))

# get sequences
insertions_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, insertions_gr)


## MSA and get consesus
# do the MSA
seq_msa <- msa::msa(insertions_seq, method = "ClustalOmega")
# seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)

# write as fasta
seq_msa@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.msa.fasta"))

# get the consensus
seq_consensus <-
  msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
  toupper(.) %>%
  Biostrings::DNAStringSet(.)
names(seq_consensus) <- "MYSERV6-int.consensus"

# save the consensus
Biostrings::writeXStringSet(seq_consensus, file = file.path(outpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.consensus.fasta"))


# # write as fasta
# seq_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("IAP.potentially_young.filtered_by_ORFs", "together", "msa.fasta", sep = ".")))
#
# # convert to seqinr format
# msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")
#
# # calculate distance matrix
# msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")
#
# # get tree and plot
# msa_tree <- ape::njs(msa_dist)
# msa_tree <- ape::ladderize(msa_tree)
#
# # save with labels
# png(filename = file.path(outpath, str_c("IAP.potentially_young.filtered_by_ORFs", "phylo_tree", "DECIPHER.msa", "labels", "png", sep = ".")), width = 2000, height = 2000, units = "px")
# plot.phylo(msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
# edgelabels(round(msa_tree$edge.length, 3), bg = "black", col = "white", font = 2)
# axisPhylo(backward = F)
# dev.off()


# ### expand ranges
# # to GRanges
# seqlevels(line1_gr) <- seqlevels(BSgenome.Cgriseus.Ensembl.CriGriPICR)
# seqlengths(line1_gr) <- seqlengths(BSgenome.Cgriseus.Ensembl.CriGriPICR)
# line1_gr_expand <- line1_gr + 10000
# line1_gr_expand <- trim(line1_gr_expand)
# 
# # get sequences
# line1_seq_expand <- getSeq(BSgenome.Cgriseus.Ensembl.CriGriPICR, line1_gr_expand)
# 
# # save as fasta
# Biostrings::writeXStringSet(x = line1_seq_expand, filepath = file.path(outpath, "LINE1.CriGriPICR.5K_to_7k.ORFs.expand_10k.fasta"))
# 
# # save
# saveRDS(line1_gr, file = file.path(outpath, str_c("LINE1_annotation.rmsk.GRanges.RDS")))
# 
# # save as .bed
# rtracklayer::export.bed(line1_gr, file.path(outpath, "LINE1_annotation.rmsk.bed"))
# 
# # save as SAF
# line1_saf <-
#   line1_gr %>%
#   as_tibble(.) %>%
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>%
#   readr::write_delim(., file.path(outpath, "LINE1_annotation.rmsk.saf"), delim = "\t")


