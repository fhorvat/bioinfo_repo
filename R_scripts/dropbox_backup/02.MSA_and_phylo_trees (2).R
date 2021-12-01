### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees/sequences")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(getwd(), "../results")

# set ensembl version
ensembl_version <- 99

# annotation with protein IDs path
annotation_proteins_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees/sequences/Piwi_Mov10l1.ensembl.proteins.csv"

# protein sequences path
protein_sequences_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees/sequences"
protein_fasta_path <- list.files(protein_sequences_path, pattern = ".*fasta$")

######################################################## READ DATA
# read annotation with protein IDs
annotation_proteins <- readr::read_csv(annotation_proteins_path)

# read protein sequences
protein_seq <- purrr::map(protein_fasta_path, function(path){
  
  # read protein sequence
  Biostrings::readAAStringSet(filepath = path)
  
}) %>% 
  set_names(basename(protein_fasta_path) %>% str_remove(., "\\.AA\\.fasta"))

######################################################## MAIN CODE
# get Mov10l and Piwil sequences in a list
protein_seq <- purrr::map(c("Mov10l1", "Piwil"), function(protein){
  
  protein_seq[str_detect(names(protein_seq), protein)] %>% 
    unname(.) %>% 
    do.call(c, .)
  
}) %>% 
  set_names(., c("Mov10l1", "Piwil"))


### phylogenetic trees  
purrr::map(c("Mov10l1", "Piwil"), function(protein){
  
  ### extract protein
  # get one protein
  one_seq <- protein_seq[[protein]]
  
  # remove human/chicken/zebrafish
  one_seq <- one_seq[!str_detect(names(one_seq), "human|chicken|zebrafish")]
  
  
  ### multiple sequence alignment and phylogeny tree
  # MSA
  protein_msa <- msa::msa(one_seq, method = "ClustalOmega", verbose = T)

  # convert to seqinr format
  protein_msa_seqinr <- msaConvert(protein_msa, type = "seqinr::alignment")
  
  # calculate distance matrix
  protein_msa_dist <- seqinr::dist.alignment(protein_msa_seqinr, "identity")
  
  # get tree and plot
  protein_msa_tree <- ape::nj(protein_msa_dist)
  protein_msa_tree <- ape::ladderize(protein_msa_tree)
  
  # save with labels
  png(filename = file.path(outpath, str_c(protein, "phylo_tree", "ClustalO.msa", "labels", "png", sep = ".")), width = 2000, height = 1000, units = "px")
  plot.phylo(protein_msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
  edgelabels(round(protein_msa_tree$edge.length, 3), bg = "black", col = "white", font = 2)
  axisPhylo(backward = F)
  dev.off()
  
  # save without labels
  png(filename = file.path(outpath, str_c(protein, "phylo_tree", "ClustalO.msa", "clean", "png", sep = ".")), width = 2000, height = 1000, units = "px")
  plot.phylo(protein_msa_tree, edge.width = 2, font = 1, show.tip.label = FALSE, cex = 1)
  dev.off()
  
  # return
  return(protein)
  
})


### Mov10l1
# get sequence
mov10l1_seq <- protein_seq[["Mov10l1"]]
mov10l1_seq <- mov10l1_seq[!str_detect(names(mov10l1_seq), "chicken|zebrafish")]

# MSA
seq_msa <- msa::msa(mov10l1_seq, method = "ClustalOmega")

# write as fasta
seq_msa@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c("Mov10l.protein", "ClustalO.msa.fasta", sep = ".")))

# # write as .pdf
# msaPrettyPrint(seq_msa,
#                output = "pdf",
#                logoColors = "chemical",
#                paperWidth = 17.5, 
#                paperHeight = 3,
#                file = file.path(outpath, str_c("Mov10l.protein", "ClustalW.msa.pdf", sep = ".")),
#                askForOverwrite = F)


### also print phylo tree
# convert to seqinr format
protein_msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")

# calculate distance matrix
protein_msa_dist <- seqinr::dist.alignment(protein_msa_seqinr, "identity")

# get tree and plot
protein_msa_tree <- ape::nj(protein_msa_dist)
# protein_msa_tree <- ape::ladderize(protein_msa_tree)
protein_msa_tree <- ape::root(protein_msa_tree, outgroup = "Mov10l1.human.ENSG00000073146.ENST00000262794")
protein_msa_tree <- ape::ladderize(protein_msa_tree)

# save with labels
png(filename = file.path(outpath, str_c("Mov10l.all_animals", "phylo_tree", "ClustalO.msa", "labels", "png", sep = ".")), width = 2000, height = 1000, units = "px")
plot.phylo(protein_msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
edgelabels(round(protein_msa_tree$edge.length, 3), bg = "black", col = "white", font = 2)
axisPhylo(backward = F)
dev.off()


# save without labels
png(filename = file.path(outpath, str_c("Mov10l.all_animals", "phylo_tree", "ClustalO.msa", "clean", "png", sep = ".")), width = 2000, height = 1000, units = "px")
plot.phylo(protein_msa_tree, edge.width = 2, font = 1, show.tip.label = FALSE, cex = 1)
dev.off()
