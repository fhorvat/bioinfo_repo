### INFO: 
### DATE: Wed Oct 10 01:31:03 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/hackaton/prepare_data/tree")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(rphast)
library(ape)
library(ggtree)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# path of nw ensembl tree
tree_path <- file.path(inpath, "species_tree.ensembl.branch_len.nw")

# species details path
species_table_path <- file.path(inpath, "../../hemoglobin.alpha_and_beta.homologs_ensembl.csv")

######################################################## READ DATA
# read Ensembl tree
tree <- rphast::read.newick.tree(tree_path)

# read species table
species_tb <- readr::read_csv(species_table_path)

######################################################## MAIN CODE
# clean species table
species_tb_tidy <- 
  species_tb %>% 
  dplyr::mutate(scientific_name = str_replace_all(scientific_name, " ", "_") %>% tolower(.))

# prune tree
pruned.tree <- rphast::prune.tree(tree, species_tb_tidy$scientific_name, all.but = TRUE)
pruned.tree <- read.tree(text = pruned.tree)

# change labels to scientific names
pruned.tree$tip.label <- species_tb_tidy$common_name[match(pruned.tree$tip.label, species_tb_tidy$scientific_name)]

# draw tree
ggtree(pruned.tree, size = 2) +
  geom_tiplab(size = 4, color = "black") +
  theme_tree2() +
  theme(legend.position = "right") +
  ggsave(filename = "species_tree.ensembl.pruned.labels.png", width = 12, height = 10) 

# # draw tree with nodes numbers

# ggtree(pruned.tree) + 
#   geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
#   geom_tiplab() + 
#   theme_tree2() +
#   theme(legend.position = "right") +
#   ggsave(filename = "species_tree.ensembl.pruned.nodes.png", width = 15, height = 10) 
