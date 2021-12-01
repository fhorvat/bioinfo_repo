### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/FLI_elements/IAP/annotation")

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
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# insertions table path
ltrs_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/LTRs/potentially_young_LTRs"
insertions_path <- file.path(ltrs_path, "potentially_young_LTRs.edit_distance.10.ORFs.csv")

######################################################## READ DATA
# read insertion table
insertions_tb <- readr::read_csv(insertions_path)

######################################################## MAIN CODE
# # tidy table, filter only IAPs, filter by ORF length (GAG = 576 aa, POL = 905 aa)
# iap_tb <- 
#   insertions_tb %>% 
#   dplyr::filter(ltr_class == "IAP") %>% 
#   dplyr::select(insertion_coordinates, strand, ltr_name, rmsk_id, longest_orf_1, longest_orf_2, insertion_class) %>% 
#   tidyr::separate(insertion_coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
#   dplyr::mutate(start = as.numeric(start), 
#                 end = as.numeric(end), 
#                 width = end - start + 1) %>% 
#   dplyr::arrange(desc(longest_orf_1 + longest_orf_2))
# 
# # save table
# readr::write_csv(iap_tb, file.path(outpath, "IAP.mm10.potentially_young.ordered_by_ORFs.20201031.csv"))

# get guys from above selected by Petr
iap_tb <- readr::read_csv(file.path(outpath, "IAP.mm10.potentially_young.ordered_by_ORFs.20201031_reduced.csv"))

# # save as SAF
# iap_tb %>% 
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
#   readr::write_delim(., file.path(outpath, "IAP.potentially_young.filtered_by_ORFs.top_49.saf"), delim = "\t")

# to GRanges
iap_gr <-
  iap_tb %>%
  GRanges(.)
names(iap_gr) <- str_c(mcols(iap_gr)$ltr_name, mcols(iap_gr)$rmsk_id, sep = ".")
names(iap_gr) <- str_replace_all(names(iap_gr), "-", "_")

# save as .bed
rtracklayer::export.bed(iap_gr, file.path(outpath, "IAP.mm10.potentially_young.ordered_by_ORFs.20201031_reduced.bed"))

# get sequences
iap_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, iap_gr)

# # save as fasta
# Biostrings::writeXStringSet(x = iap_seq, filepath = file.path(outpath, "IAP.mm10.potentially_young.ordered_by_ORFs.20201031_reduced.fasta"))


### get consensus by repeat names
# split by repeat names
iap_seq_list <- split(iap_seq, iap_gr$ltr_name)

## MSA and get consesus
iap_seq_consensus <- purrr::map(names(iap_seq_list), function(seq_name){
  
  # do the MSA
  seq_msa <- msa::msa(iap_seq_list[[seq_name]], method = "ClustalOmega")
  # seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)
  
  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("IAP.mm10.potentially_young.ordered_by_ORFs.20201031_reduced", seq_name, "msa.fasta", sep = ".")))
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # return
  return(seq_consensus)
  
})

# join together
iap_seq_consensus <- do.call(c, iap_seq_consensus)
names(iap_seq_consensus) <- names(iap_seq_list)

# save the consensus
Biostrings::writeXStringSet(iap_seq_consensus, file = file.path(outpath, str_c("IAP.mm10.potentially_young.ordered_by_ORFs.20201031_reduced", "consensus.fasta", sep = ".")))


# ## MSA together
# # MSA using DECIPHER
# seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq, iterations = 100, refinements = 100, processors = 1, verbose = T)
# 
# # msa
# seq_msa_2 <- msa::msa(line1_seq, method )
# 
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

