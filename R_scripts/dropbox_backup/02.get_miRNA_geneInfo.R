### INFO: gets miRNA annotation
### DATE: Fri Jun 04 16:23:06 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/annotation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

######################################################## READ DATA
### get miRNA annotations
# set genome list
genome_list <- c("cow.bosTau9", "ghamster.mesAur1", "human.hg38", "mouse.mm10", "pig.susScr11")

# map through genomes 
purrr::map(genome_list, function(genome){
  
  # link path
  gtf_path <- list.files(inpath, genome, full.names = T)
  
  # get full link
  gtf_path_full <- Sys.readlink(gtf_path)
  
  # read gtf, get gene info
  gtf_tb <- 
    rtracklayer::import.gff(gtf_path_full) %>% 
    as_tibble(.) %>% 
    dplyr::select(gene_id = ID, seqnames, start, end, strand,
                  gene_name = Name, gene_biotype = type) %>% 
    dplyr::mutate(gene_description = NA)
  
  # save
  geneInfo_path <- str_replace(gtf_path_full, "\\.gff.*$", ".geneInfo.csv")
  if(geneInfo_path != gtf_path_full){
    readr::write_csv(gtf_tb, geneInfo_path)
  }else{
    stop("You're trying to overwrite your .gtf file, please stop")
  }
  
  
})
