### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/map_mouse_chromosomes_to_Siomi")

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
library(BSgenome.Maur.UCSC.MesAur1)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped clusters path
bam_path <- file.path(inpath, "Siomi_to_mm10.bam")

######################################################## READ DATA
# get bam
bam <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T)

######################################################## MAIN CODE
# get longest alignment for each chromosome
bam_gr <- 
  bam %>% 
  grglist(.) %>% 
  range(.) %>% 
  unlist(.)

# get names as coordinates
mcols(bam_gr)$mm10.chr <- names(bam_gr)

# get in table
bam_tb <- 
  bam_gr %>% 
  as_tibble(.) %>% 
  dplyr::group_by(seqnames) %>% 
  dplyr::filter(width == max(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(seqnames)

