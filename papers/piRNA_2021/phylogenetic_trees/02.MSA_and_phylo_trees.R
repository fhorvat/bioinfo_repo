### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees")

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

library(biomaRt)
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
inpath <- getwd()

# set outpath
outpath <- getwd()

# annotation with protein IDs path
annotation_proteins_path <- file.path(inpath, "animal_list.protein_ids.20201014.csv")

# protein sequences path
protein_sequences_path <- file.path(inpath, "Mov10l1_sequences")
protein_fasta_path <- list.files(protein_sequences_path, pattern = ".*fasta$", full.names = T)
protein_fasta_path <- protein_fasta_path[!str_detect(protein_fasta_path, "tunicate")]

######################################################## READ DATA
# read annotation with protein IDs
annotation_proteins <- readr::read_csv(annotation_proteins_path)

# read protein sequences
protein_seq <- Biostrings::readAAStringSet(filepath = protein_fasta_path)
names(protein_seq) <- basename(protein_fasta_path) %>% str_remove(., "\\..*")

######################################################## MAIN CODE
### multiple sequence alignment and phylogeny tree
# MSA with msa package
protein_msa_ClustalO <- msa::msa(protein_seq, method = "ClustalOmega", verbose = T)

# MSA with DECIPHER package
protein_msa_decipher <- 
  DECIPHER::AlignSeqs(myXStringSet = protein_seq, iterations = 100, refinements = 100, processors = 1, verbose = T) %>% 
  AAMultipleAlignment(.)

# create an list
protein_msa_list <- list("ClustalO" = protein_msa_ClustalO, "DECIPHER" = protein_msa_decipher)

### phylogeny tree and plot
# do it for msa and DECIPHER
purrr::map(names(protein_msa_list), function(software){
  
  software <- "ClustalO"
  
  # get MSA from one software
  msa_protein <- protein_msa_list[[software]]
  
  # write as fasta
  msa_protein@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("Mov10l1.protein", software, "fasta", sep = ".")))
  
  # # write as .pdf
  # msaPrettyPrint(seq_msa,
  #                output = "pdf",
  #                logoColors = "chemical",
  #                paperWidth = 17.5, 
  #                paperHeight = 3,
  #                file = file.path(outpath, str_c("Mov10l.protein", "ClustalW.msa.pdf", sep = ".")),
  #                askForOverwrite = F)
  
  # convert to seqinr format
  msa_seqinr <- msaConvert(msa_protein, type = "seqinr::alignment")
  
  # calculate identity distance matrix 
  # identity = ratio of the number of matching residues to the total length of the alignment) 
  # for example, if identity between 2 sequences is 80 it's the squared root of (1.0 - 0.8) = 0.4472136
  msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")
  
  # # calculate sequence difference
  msa_dist <- msa_dist^2

  # get tree and plot
  msa_tree <- ape::nj(msa_dist)
  # msa_tree <- ape::root(msa_tree, outgroup = "shark")
  # msa_tree <- ape::ladderize(msa_tree)
  
  # plot using ggtree
  ggtree(msa_tree) +
    geom_tiplab(geom = "label", hjust = .5) +
    geom_label(aes(x = branch, label = round(branch.length, 3)), fill = "white") +
    scale_y_reverse() + 
    theme_tree2() +
    ggsave(file.path(outpath, str_c("Mov10l1.protein", "phylo_tree", software, "labels", "ggtree", "png", sep = ".")), width = 20, height = 10)
  
  # # save with labels
  # png(filename = file.path(outpath, str_c("Mov10l1.protein", "phylo_tree", software, "labels", "png", sep = ".")), width = 2000, height = 2000, units = "px")
  # plot.phylo(msa_tree, edge.width = 4, font = 4, show.tip.label = TRUE, cex = 2, use.edge.length = TRUE)
  # edgelabels(round(msa_tree$edge.length, 3), bg = "black", col = "white", font = 4, cex = 2)
  # # axisPhylo(backward = F)
  # dev.off()
  
  # # save without labels
  # png(filename = file.path(outpath, str_c("Mov10l.protein", "phylo_tree", software, "png", sep = ".")), width = 1500, height = 1200, units = "px")
  # plot.phylo(msa_tree, edge.width = 6, font = 1, show.tip.label = FALSE, cex = 1)
  # dev.off()
  
  # return
  return(software)
  
})


