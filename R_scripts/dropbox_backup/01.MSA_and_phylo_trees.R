### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Mollusca_project/orhtologs")

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

library(taxize)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# protein table path
protein_path <- file.path(inpath, "Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences.xlsx")

######################################################## READ DATA
# read protein table
protein_tb <- read_xlsx(protein_path)

# get non-redundant list of animals
protein_tb %>% 
  dplyr::select(organism_name, phylum, subphylum, class) %>% 
  dplyr::distinct(.) %>% 
  dplyr::mutate(phylum = factor(phylum, levels = c("Placozoa", "Porifera", "Cnidaria",
                                                   "Echinodermata", "Hemichordata", "Chordata", 
                                                   "Mollusca", "Annelida", "Brachiopoda", "Platyhelminthes", 
                                                   "Priapulida", "Nematoda", "Arthropoda"))) %>% 
  dplyr::arrange(phylum, subphylum, class, organism_name) %>% 
  openxlsx::write.xlsx(., file = file.path(outpath, "Dicer1.Metazoa.OrthoDB.20211117.animal_list.xlsx"))

######################################################## MAIN CODE
# filter table
protein_tb <-
  protein_tb %>%
  dplyr::group_by(organism_name) %>%
  dplyr::slice_max(order_by = nchar(protein_seq)) %>%
  dplyr::ungroup(.)
#   dplyr::group_by(phylum) %>% 
#   dplyr::slice_sample(n = 10) %>% 
#   dplyr::ungroup(.)

# get protein sequences as AAStringSet
protein_seq <- Biostrings::AAStringSet(protein_tb$protein_seq)
names(protein_seq) <- protein_tb$int_prot_id

# set output name
output_name <- "Dicer1.Metazoa.OrthoDB.20211117"

# ### multiple sequence alignment and phylogeny tree
# # MSA with msa package
# msa_protein <- msa::msa(protein_seq, method = "ClustalOmega", verbose = T)

# MSA with DECIPHER package
msa_protein <- 
  DECIPHER::AlignSeqs(myXStringSet = protein_seq, iterations = 100, refinements = 100, processors = 1, verbose = T) %>% 
  AAMultipleAlignment(.)

# set software name
software <- "DECIPHER"

# write as fasta
msa_protein@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c(output_name, "MSA", software, "fasta", sep = ".")))

# convert to seqinr format
msa_seqinr <- msaConvert(msa_protein, type = "seqinr::alignment")

# calculate identity distance matrix 
# identity = ratio of the number of matching residues to the total length of the alignment) 
# for example, if identity between 2 sequences is 80 it's the squared root of (1.0 - 0.8) = 0.4472136
msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")

# # calculate sequence difference
msa_dist <- msa_dist^2

# get tree and plot
msa_tree <- ape::njs(msa_dist)
# msa_tree <- ape::root(msa_tree, outgroup = "shark")
# msa_tree <- ape::ladderize(msa_tree)

# get tree nodes and names in table, join with annotation
tree_tb <- 
  tibble(node = nodeid(msa_tree, names(protein_seq)),
         int_prot_id = names(protein_seq)) %>% 
  dplyr::left_join(., protein_tb, by = "int_prot_id") %>% 
  dplyr::mutate(phylum = factor(phylum, levels = c("Placozoa", "Porifera", "Cnidaria",
                                                   "Echinodermata", "Hemichordata", "Chordata", 
                                                   "Mollusca", "Annelida", "Brachiopoda", "Platyhelminthes", 
                                                   "Priapulida", "Nematoda", "Arthropoda")))
  # dplyr::mutate(tree_color = gg_color_hue(13)[as.numeric(phylum)])

# join with tree
msa_tree_annotated <- full_join(msa_tree, tree_tb, by = "node")

# plot using ggtree
tree_plot <-
  ggtree(msa_tree_annotated, 
         aes(color = phylum),
         layout = "fan", 
         ladderize = TRUE, continuous = FALSE, size = 1) +
  scale_color_discrete() +
  # geom_tiplab(aes(angle = angle,
  #                 label = organism_name), 
  #             color = "black") +
  # geom_tippoint(aes(color = phylum), size = 5) +
  theme_tree2() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()) + 
  ggsave(file.path(outpath, str_c(output_name, "phylo_tree", software, "ggtree", "png", sep = ".")), 
         width = 15, height = 15)

