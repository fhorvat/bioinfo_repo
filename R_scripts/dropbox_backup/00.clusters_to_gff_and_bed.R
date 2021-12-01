### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/small_RNAseq")

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

library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# cluster table path
cluster_tb_list <- list.files(inpath, "piRNA_clusters.*\\.20210513\\.csv", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
### loop through tables
purrr::map(cluster_tb_list, function(cluster_tb_path){
  
  # read cluster table
  cluster_tb <- readr::read_csv(cluster_tb_path)
  
  # get cluster coordinates
  cluster_gr <- 
    cluster_tb %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
    dplyr::mutate(strand = "*") %>% 
    GRanges(.)
  
  # set names
  names(cluster_gr) <- str_c(seqnames(cluster_gr), ":", start(cluster_gr), "-", end(cluster_gr))
  
  # export as .gtf
  export.gff3(object = cluster_gr, con = str_replace(cluster_tb_path, "\\.csv$", ".gff"))
  
  # export as .bed
  export.bed(object = cluster_gr, con = str_replace(cluster_tb_path, "\\.csv$", ".bed"))
  
})

