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

####################################################### READ DATA
# # set animals which to include in tree
# animals <- c("Mouse",
#              "Rat",
#              "Kangaroo_rat",
#              "Naked_mole_rat",
#              "Guinea_pig",
#              "Squirrel",
#              "Rabbit",
#              "Cow",
#              "Human"
# )
# 
# # read UCSC tree
# tree <- read.newick.tree("mm10.60way.commonNames.nh")
# 
# # prune tree to selected animals
# pruned.tree <- prune.tree(tree, animals, all.but = TRUE)

######################################################## MAIN CODE
# put human outside, convert tree, add labels
pruned.tree <- "(((((((Mouse:0.084509,Rat:0.091589):0.197773,Kangaroo_rat:0.211609):0.022992,(Naked_mole_rat:0.1,Guinea_pig:0.125629):0.1):0.01015,Squirrel:0.148468):0.025746,Rabbit:0.21569):0.015313,Cow:0.18908):0.020593,Human:0.07);"
pruned.tree <- read.tree(text = pruned.tree)
pruned.tree$tip.label <- c("Mouse", "Rat",
                           "Kangaroo rat", "Naked mole rat", "Guinea pig", 
                           "Squirrel", "Rabbit",        
                           "Cow", "Human")

## set intron and insert
# intron 1 - ORR1D2 (mm10: chr15:73178442-73178690)
# intron 1 - ORR1B1 (mm10: chr15:73159826-73160209)
# intron 3 - MTD (mm10: chr15:73132252-73132555)
# intron 3 - MTD (mm10:  chr15:73134556-73134811)
# intron 3 - ORR1B1 (mm10:  chr15:73135071-73135351)
# intron 3 - MTA_Mm (mm10:  chr15:73135562-73135637)
introns <- c("intron 1", "intron 1", "intron 3", "intron 3", "intron 3", "intron 3")
inserts <- c("ORR1D2", "ORR1B1", "MTD", "MTD", "ORR1B1", "MTA_Mm")
insert_poss <- c("mm10: chr15:73178442-73178690", "mm10: chr15:73159826-73160209", "mm10: chr15:73132252-73132555", 
                 "mm10: chr15:73134556-73134811", "mm10: chr15:73135071-73135351", "mm10: chr15:73135562-73135637")
which_animals_list <- list(c("Mouse", "Rat"), 
                           c("Mouse", "Rat"), 
                           c("Mouse", "Rat"), 
                           c("Mouse", "Rat"), 
                           c("Mouse", "Rat"), 
                           c("Mouse"))


### group and draw tree in loop
for(n in 1:length(introns)){
  
  intron <- introns[n]
  insert <- inserts[n]
  insert_pos <- insert_poss[n]
  which_animals <- which_animals_list[[n]]
  
  # group tree
  pruned.tree <- groupOTU(object = pruned.tree, 
                          focus = which_animals)
  
  # relevel
  attr(pruned.tree, "group") <- relevel(attr(pruned.tree, "group"), ref = 2)
  
  # draw tree, rotate
  draw_tree <- 
    ggtree(pruned.tree, aes(color = group)) +
    geom_tiplab(size = 3, color = "black")
  draw_tree <- rotate(draw_tree, 11)
  draw_tree <- rotate(draw_tree, 10)
  draw_tree <- rotate(draw_tree, 16)

  ## label tree and save
  draw_tree + 
    # geom_text2(aes(subset = !isTip, label = node), hjust = - .3) +
    scale_color_manual(labels = c(insert, str_c("no ", insert)), values = c("red", "gray60")) + 
    theme_tree2() +
    labs(title = str_c(intron, ", ", insert, ", ", insert_pos), color = "Insert") +
    theme(legend.position = "right") +
    ggsave(filename = str_c("phyloTree.noHamsters.", 
                            str_replace_all(intron, " ", "_"), ".", 
                            insert, ".", 
                            str_replace_all(insert_pos, "mm10: ", ""), 
                            ".pdf"), width = 12, height = 5) 
  
}
