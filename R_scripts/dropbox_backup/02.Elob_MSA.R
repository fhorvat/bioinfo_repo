### INFO: 
### DATE: Tue Nov 26 14:58:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Elob_MSA")
# setwd("C:/Users/fhorvat/Dropbox/Bioinfo/Svoboda/lnc1_paper/Elob_MSA")

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
library(Biostrings)
library(msa)
library(seqinr)
library(ape)
library(ggtree)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "sequences")

# set outpath
outpath <- getwd()

# protein sequences path
protein_fasta_path <- list.files(inpath, pattern = "Elobl*.protein\\..*\\.fa", full.names = T)

######################################################## READ DATA
# read protein fasta
protein_seq <- Biostrings::readAAStringSet(protein_fasta_path)
names(protein_seq) <- str_remove(names(protein_seq), " peptide: .*")

######################################################## MAIN CODE
### multiple sequence alignment and phylogeny tree
# MSA
seq_msa <- msa::msa(protein_seq, method = "ClustalOmega")

# write as fasta
seq_msa@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c("Elob_Elobl.protein", "ClustalO.msa.fasta", sep = ".")))

# write as .pdf
msaPrettyPrint(seq_msa,
               output = "pdf",
               logoColors = "chemical",
               paperWidth = 17.5, paperHeight = 3,
               file = file.path(outpath, str_c("Elob_Elobl.protein", "ClustalO.msa.pdf", sep = ".")),
               askForOverwrite = F)
# tools::texi2pdf(file.path(outpath, str_c("Elob_Elobl.protein", "ClustalO.msa.tex", sep = ".")))

### tree
# convert to seqinr format
seq_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")

# calculate distance matrix
seq_msa_dist <- seqinr::dist.alignment(seq_seqinr, "identity")

# get tree
seq_msa_tree <- ape::nj(seq_msa_dist)
# seq_msa_tree <- ape::ladderize(seq_msa_tree)
seq_msa_tree$tip.label <- str_remove(seq_msa_tree$tip.label, " .*")

# plot
ggtree(seq_msa_tree) +
  geom_tiplab(geom = "label", hjust = .5) +
  # xlim(NA, 8) +
  theme(legend.position = "right") +
  scale_color_viridis_c() +
  ggsave(file.path(outpath, "Elob_Elobl.protein.MSA_tree.png"), width = 10, height = 5)
    

