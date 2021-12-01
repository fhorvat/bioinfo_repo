#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 24.03.2017. 
### AUTHOR: Filip Horvat

################################################################################### SCRIPT PARAMS
rm(list = ls()); gc()
options(bitmapType = "cairo")

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/paired_end/mouse_strains/Data/Mapped/STAR_mm10_merged")

################################################################################### LIBRARIES
# data shaping
library(magrittr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(scales)

################################################################################### SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

################################################################################### FUNCTIONS

################################################################################### PATH VARIABLES
inpath <- getwd()
outpath <- getwd()

################################################################################### TABLES
# read histogram table from bbmerge.sh output
hist_df <- read_delim(file = file.path(inpath, "s_SG_BC7_Cast2.PE.merged.hist.txt"), delim = "\t", col_names = T, comment = "#")

################################################################################### MAIN CODE
# plot histogram
ggplot(data = hist_df) +
  geom_col(aes(x = InsertSize, y = Count), color = "grey30", fill = "grey30") +
  scale_y_continuous(labels = comma) +
  xlab("insert size") +
  ylab("read count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15), 
        axis.text.y = element_text(hjust = 1, size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(filename = "s_SG_BC7_Cast2.PE.merged.hist.png", width = 15, height = 10)

