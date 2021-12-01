### INFO: reads .out files from repeatMasker
### DATE: 13. 02. 2018.
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/trees")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(rphast)
library(ape)
library(ggtree)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read data.frame with data about animals
animals_sci <- 
  readr::read_csv(file = "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/Ago2_mammals.csv") %>% 
  dplyr::mutate(sci_name_tree = scientific_name %>% tolower(.) %>% stringr::str_replace_all(., " ", "_"),
                sci_name_tree = replace(sci_name_tree, sci_name_tree == "mus_spretus", "mus_spretus_spreteij"), 
                sci_name_tree = replace(sci_name_tree, sci_name_tree == "heterocephalus_glaber", "heterocephalus_glaber_male"), 
                sci_name_tree = replace(sci_name_tree, sci_name_tree == "cricetulus_griseus", "cricetulus_griseus_crigri"))

# # read Ensembl tree
# tree <- read.newick.tree("species_tree.67mammals.branch_len.nw")
# pruned.tree <- prune.tree(tree, animals_sci$sci_name_tree, all.but = TRUE)
# pruned.tree <- read.tree(text = pruned.tree)

# read tree with added Apodemus sylvaticus
pruned.tree <- read.tree("Ago2.animals.plus.asylvaticus.nw")

######################################################## MAIN CODE
# change labels to scientific names
pruned.tree$tip.label <- animals_sci$sci_name[match(pruned.tree$tip.label, animals_sci$sci_name_tree)]

## set intron and insert
# intron 1 - ORR1D2 (mm10: chr15:73178442-73178690) + 
# intron 3 - MTD (mm10:  chr15:73134556-73134811) + 
# intron 3 - MTA_Mm (mm10:  chr15:73135562-73135637) +
introns <- c("intron 1", "intron 3", "intron 3")
inserts <- c("ORR1D2", "MTD", "MTA_Mm")
insert_poss <- c("mm10: chr15:73178442-73178690", "mm10: chr15:73134556-73134811", "mm10: chr15:73135562-73135637")
which_animals_list <- list(c("Mus musculus", "Mus spretus", "Apodemus sylvaticus", "Rattus norvegicus", 
                             "Cricetulus griseus", "Mesocricetus auratus", "Microtus ochrogaster", "Peromyscus maniculatus bairdii",
                             "Jaculus jaculus"), 
                           c("Mus musculus", "Mus spretus", "Apodemus sylvaticus", "Rattus norvegicus",
                             "Cricetulus griseus", "Mesocricetus auratus", "Microtus ochrogaster", "Peromyscus maniculatus bairdii"), 
                           c("Mus musculus", "Mus spretus"))

### group and draw tree in loop
for(n in 1:length(introns)){
  
  intron <- introns[n]
  insert <- inserts[n]
  insert_pos <- insert_poss[n]
  which_animals <- which_animals_list[[n]]
  
  # group tree
  pruned.tree <- groupOTU(object = pruned.tree, focus = which_animals)
  
  # relevel
  attr(pruned.tree, "group") <- relevel(attr(pruned.tree, "group"), ref = 2)
  
  # draw tree, rotate
  draw_tree <- 
    ggtree(pruned.tree, aes(color = group)) +
    geom_tiplab(size = 2.5, color = "black")
  
  # draw_tree <- rotate(draw_tree, 13)
  # draw_tree <- rotate(draw_tree, 12)
  # draw_tree <- rotate(draw_tree, 19)
  # draw_tree <- rotate(draw_tree, 18)
  
  ## label tree and save
  draw_tree + 
    # geom_text2(aes(subset = !isTip, label = node), hjust = - .3) +
    geom_cladelabel(node = 3, label = "", align = T, offset = 0.03, fontsize = 2) +
    geom_cladelabel(node = 3, label = "Sciuridae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 14, label = "Heterocephalidae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 15, label = "Caviidae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 4, label = "Heteromyidae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 5, label = "Dipodidae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 28, label = "Cricetidae", align = T, offset = 0.017, fontsize = 3) +
    geom_cladelabel(node = 25, label = "Muridae", align = T, offset = 0.017, fontsize = 3) +
    scale_color_manual(labels = c(insert, str_c("no ", insert)), values = c("red", "gray60")) + 
    theme_tree2() +
    labs(title = str_c(intron, ", ", insert, ", ", insert_pos), color = "Insert") +
    theme(legend.position = "right") +
    ggsave(filename = str_c("phyloTree.", 
                            str_replace_all(intron, " ", "_"), ".", 
                            insert, ".", 
                            str_replace_all(insert_pos, "mm10: ", ""), 
                            ".pdf"), width = 12, height = 5) 
  
}
