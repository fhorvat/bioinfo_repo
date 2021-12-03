### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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
library(DECIPHER)
library(seqinr)
library(ape)
library(ggtree)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# protein table path
protein_path <- file.path(inpath, "..", "Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences_PS_CDD.xlsx")

# tree path
tree_path <- list.files(inpath, pattern = "RAxML_bestTree.*supertree.raxml", full.names = T)

######################################################## READ DATA
# read protein table
protein_tb <- read_xlsx(protein_path)

# read tree
msa_tree <- ape::read.tree(tree_path)

######################################################## MAIN CODE
# add factor levels
protein_tb %<>%
  dplyr::mutate(phylum = factor(phylum, levels = c("Placozoa", "Porifera", "Cnidaria",
                                                   "Echinodermata", "Hemichordata", "Chordata",
                                                   "Mollusca", "Annelida", "Brachiopoda", "Platyhelminthes",
                                                   "Priapulida", "Nematoda", "Arthropoda"))) %>%
  dplyr::arrange(phylum, subphylum, class, organism_name, desc(protein_length)) %>% 
  dplyr::mutate(class_name = ifelse(is.na(class), subphylum, class), 
                class_name = ifelse(is.na(class_name), phylum, class_name), 
                class_name = factor(class_name, levels = unique(class_name))) %>% 
  dplyr::mutate(protein_name = str_replace_all(organism_name, " ", "_") %>% 
                  str_c(., pub_gene_id, sep = "_") %>% 
                  str_replace_all(., ";", "_"))

# get tree name
tree_name <- 
  tree_path %>% 
  basename(.)
# str_remove(., "\\.tree$")

# get tree nodes and names in table, join with annotation
tree_tb <- 
  tibble(node = nodeid(msa_tree, protein_tb$protein_name),
         protein_name = protein_tb$protein_name) %>% 
  dplyr::left_join(., protein_tb, by = "protein_name")

# join with tree
msa_tree_annotated <- full_join(msa_tree, tree_tb, by = "node")

# plot using ggtree
tree_plot <-
  ggtree(msa_tree, 
         # layout = "fan", 
         ladderize = TRUE, 
         continuous = FALSE,
         size = 1) +
  scale_color_discrete() +
  theme_tree2() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), 
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 12)))

# get xlim
x_lim <- max(ggplot_build(tree_plot)$data[[3]]$x)

# extend x-limit 
tree_plot <- 
  tree_plot + 
  xlim(0, x_lim + 1)

# save
ggsave(plot = tree_plot, 
       filename = file.path(outpath, str_c(tree_name, "ggtree", "pdf", sep = ".")), 
       width = 15, height = 15)



