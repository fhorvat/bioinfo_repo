### INFO: 
### DATE: Fri Jul 23 11:19:33 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/mouse_testis.Papd7.small_RNAseq.Apr_2021/Analysis/read_class")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# all grid files
grid_path <- list.files(inpath, "grid_[0-9]+\\.RDS")

######################################################## READ DATA
# read grid files
grid_list <- purrr::map(grid_path, readRDS)

######################################################## MAIN CODE
# set experiment name
experiment_name <- "mouse_testis.Papd7.small_RNAseq.Apr_2021"

# # plot on one page
# grid <- 
#   plot_grid(grid_list[[1]],grid_list[[2]], 
#             nrow = 2) + 
#   ggsave(filename = file.path(outpath, "all.pdf"), 
#          width = 20, 
#          height = 24)

# plot on one page
pdf(str_c(experiment_name, "reads_overview.pdf", sep = "."), width = 25, height = 20)
invisible(lapply(grid_list, print))
dev.off()
