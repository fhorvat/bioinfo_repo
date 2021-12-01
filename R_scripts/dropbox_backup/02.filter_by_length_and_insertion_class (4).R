### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation")

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
library(msa)

library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.FIXED"

# repeatMasker
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# LINE1s rmsk path
# line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.annotated_exons.csv")
line1_path <- file.path(inpath, "LINE1.4000nt_plus.ORFs.csv")

######################################################## READ DATA
# read LINE1 table
line1_tb <- readr::read_csv(line1_path)

# read repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# filter table - ORF1 at least 370 amino-acids long and ORF2 should be at least 1200 amino acids long
line1_out <- 
  line1_tb %>% 
  # dplyr::filter(width >= 5000, width <= 7000) %>% 
  # dplyr::filter(insertion_class %in% c("whole", "within")) %>% 
  dplyr::filter(longest_orf_1 >= 1200*3, 
                longest_orf_2 >= 370*3) %>% 
  dplyr::arrange(-longest_orf_1) %>% 
  dplyr::select(-c(repClass, repFamily)) %>% 
  dplyr::mutate(repName = repName %>% 
                  str_split(., "/") %>% 
                  purrr::map(., function(name) unique(name) %>% str_c(., collapse = "/")) %>% 
                  unlist(.))

# save table
line1_out %T>%
  readr::write_csv(., file.path(outpath, "LINE1.Siomi.with_both_ORFs.csv"))

# to GRanges
line1_gr <-
  line1_out %>%
  GRanges(.)
names(line1_gr) <- str_c(mcols(line1_gr)$repName, mcols(line1_gr)$rmsk_id, sep = ".")
names(line1_gr) <- str_replace_all(names(line1_gr), "-", "_")

# get sequences
line1_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, line1_gr)

# save as fasta
Biostrings::writeXStringSet(x = line1_seq, filepath = file.path(outpath, "LINE1.Siomi.with_both_ORFs.fasta"))


### get consensus by repeat names
# split by repeat names
line1_seq_list <- split(line1_seq, line1_gr$repName)

## MSA and get consesus
line1_seq_consensus <- purrr::map(names(line1_seq_list), function(seq_name){

  # do the MSA
  seq_msa <- msa::msa(line1_seq_list[[seq_name]], method = "ClustalOmega")

  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("LINE1.Siomi.with_both_ORFs", seq_name, "msa.fasta", sep = ".")))

  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)

  # return
  return(seq_consensus)

})

# join together
line1_seq_consensus <- do.call(c, line1_seq_consensus)
names(line1_seq_consensus) <- names(line1_seq_list)

# save the consensus
Biostrings::writeXStringSet(line1_seq_consensus, file = file.path(outpath, str_c("LINE1.Siomi.with_both_ORFs", "consensus.fasta", sep = ".")))


# ## MSA together
# # msa
# seq_msa <- msa::msa(line1_seq, method = "ClustalOmega")
# 
# # write as fasta
# seq_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("LINE1.Siomi.with_both_ORFs", "together", "msa.fasta", sep = ".")))
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
# png(filename = file.path(outpath, str_c("LINE1.Siomi.with_both_ORFs", "phylo_tree", "ClustalO.msa", "labels", "png", sep = ".")), width = 2000, height = 2000, units = "px")
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

