### INFO: 
### DATE: Sat Sep 08 21:21:46 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# ESC clusters path
clust_ESC_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/ES_DcrTrans_2012/clusters.ESC.DcrO.DcrS.union.rpm.3.width.50.csv"

# 3T3 clusters path
clust_3T3_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/T3T_DcrTrans_2011/clusters.3T3.DicerO.DicerS.union.rpm.3.width.50.csv"

# annotations path 
annot_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/l.annot.RDS"

######################################################## READ DATA
# read ESC clusters
clust_ESC <- readr::read_csv(clust_ESC_path)

# read 3T3 clusters
clust_3T3 <- readr::read_csv(clust_3T3_path)

# read annotation
l.annot <- readRDS(annot_path)

######################################################## MAIN CODE
# get Mos genomic coordinates
mos_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000078365")]

## get GRanges of clusters, subset by Mos
# ESC
clust_ESC_gr <- 
  clust_ESC %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GenomicRanges::GRanges(.) %>% 
  IRanges::subsetByOverlaps(., mos_gr) %>% 
  as.data.frame(.) %>% 
  as.tibble(.)

# 3T3
clust_3T3_gr <- 
  clust_3T3 %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GenomicRanges::GRanges(.) %>% 
  IRanges::subsetByOverlaps(., mos_gr) %>% 
  as.data.frame(.) %>% 
  as.tibble(.)

