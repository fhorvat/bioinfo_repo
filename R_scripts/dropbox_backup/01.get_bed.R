### INFO: 
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Data/Mapped/STAR_DicerLocus/STAR_index")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# annotation
annot_gr <- 
  tibble(seqnames = c("Dicer1.exon2_exon3", "Dicer1.exon2_exon3", "Dicer1.exon2_HAtag_exon7", "Dicer1.exon2_HAtag_exon7", "Dicer1.exon2_HAtag_exon7"),
         start = c(1, 188, 1, 80, 107), 
         end = c(188, 351, 80, 107, 275), 
         pos = c("exon2", "exon3", "exon2", "HAtag", "exon7")) %>% 
  GRanges(.) 

names(annot_gr) <- annot_gr$pos
# mcols(annot_gr)$pos <- NULL

rtracklayer::export.bed(annot_gr, file.path(outpath, "Dicer1.locus.bed"))


