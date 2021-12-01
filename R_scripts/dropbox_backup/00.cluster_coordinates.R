### INFO: 
### DATE: Sat Sep 08 15:03:19 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Analysis/expression/smallRNA_clusters/diffExp")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(data.table)
library(GenomicRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# cluster expression path
cluster_path <- file.path(inpath, "clusters.DicerX_embryos.WT.KO.union.rpm.3.width.50.mean_RPM.detailed.csv")

######################################################## READ DATA
# read cluster expression
cluster_df <- readr::read_csv(cluster_path)

######################################################## MAIN CODE
# get clusters from expression table
clust_all <- 
  cluster_df %>% 
  dplyr::select(coordinates, strand, class, hit) %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(gene_id = str_c("cl.", 1:n())) %>% 
  GenomicRanges::GRanges(.)

# write as .gtf
rtracklayer::export.gff3(clust_all, file.path(outpath, "clusters.DicerX_embryos.WT.KO.union.rpm.3.width.50.gff3"))

