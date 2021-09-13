### INFO: 
### DATE: Tue Aug 03 14:21:57 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/nontemplated_additions")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list of datasets
datasets_list <- c("mouse.mm10.GarciaLopez_2015_RNA_GSE59254", "cow.bosTau9", "pig.susScr11")
datasets_list <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/expression",
                           datasets_list)

# expression paths
fpm_mean_path <- list.files(datasets_list, ".*\\.miRNA\\.FPM_mean\\.csv", full.names = T)

######################################################## READ DATA
# read expression
fpm_mean_list <- purrr::map(fpm_mean_path, function(path){
  
  readr::read_csv(path) %>% 
    dplyr::select_at(vars("gene_id", "gene_name", "coordinates", "strand", matches("GV|oocytes|oocyte_large"))) %>% 
    dplyr::rename_at(vars(matches("GV|oocytes|oocyte_large")), ~("FPM_oocyte"))
  
}) %>% 
  set_names(., str_extract(fpm_mean_path, "mouse.mm10.GarciaLopez_2015_RNA_GSE59254|cow.bosTau9|pig.susScr11"))

######################################################## MAIN CODE
# get the to 10 most abundant miRNAs
mirna_top <- purrr::map(names(fpm_mean_list), function(animal){
  
  fpm_mean_list[[animal]] %>% 
    dplyr::top_n(n = 10, wt = FPM_oocyte) %>% 
    dplyr::arrange(-FPM_oocyte)
  
}) %>% 
  set_names(., names(fpm_mean_list))

# get coordinates, save as bed
purrr::map(names(mirna_top), function(animal){
  
  # save info
  mirna_top[[animal]] %T>% 
    readr::write_csv(., file = file.path(outpath, str_c(animal, "top10_miRNA", "geneInfo.csv", sep = ".")))
  
  # get GRanges
  mirna_top_coords <- 
    mirna_top[[animal]] %>% 
    dplyr::select(coordinates, strand, gene_id) %>% 
    tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
    GRanges(.)
  
  # set names
  names(mirna_top_coords) <- mcols(mirna_top_coords)$gene_id
  mcols(mirna_top_coords) <- NULL
  
  # save
  rtracklayer::export.bed(mirna_top_coords, file.path(outpath, str_c(animal, "top10_miRNA", "bed", sep = ".")))
  
  # return
  return(mirna_top_coords)
  
}) %>% 
  set_names(., names(mirna_top))


