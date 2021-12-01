### INFO: gets 3'UTR sequences, longest for each gene
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/miR-205_pig/oocyte_seed_expression")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# FPM table path
fpm_path <- list.files(inpath, ".*\\.FPM_mean\\.csv", full.names = T, recursive = T)

# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# set ensembl version
ensembl_version <- 99

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# repeatMasker path
rmsk_tb_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

######################################################## READ DATA
# read FPM tables
fpm_tb_list <- purrr::map(fpm_path, function(path){
  
  # read
  readr::read_csv(path)
  
}) %>% 
  set_names(., str_extract(fpm_path, "7mer|8mer"))

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read repeatMasker table
rmsk_tb <- readr::read_delim(rmsk_tb_path, delim = "\t")

######################################################## MAIN CODE
### clean files
# get repeatMasker GRanges
rmsk_gr <- GRanges(rmsk_tb)

# loop for 7mer and 8mer
purrr::map(names(fpm_tb_list), function(name){
  
  # get table, filter top 100, get coordinates, save
  fpm_tb <- 
    fpm_tb_list[[name]] %>%
    dplyr::rename(FPM = GV) %>% 
    dplyr::top_n(x = ., n = 100, wt = FPM) %>% 
    dplyr::arrange(-FPM) %>% 
    dplyr::mutate(coordinates = str_replace_all(gene_id, ":|-", " ")) %>% 
    dplyr::select(coordinates, FPM)
  
  # get GRanges
  fpm_gr <- 
    fpm_tb %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
    GRanges(.)
  
  ### overlap with ensembl
  overlap_ensembl <- findOverlaps(fpm_gr, exons_gr, ignore.strand = T)
  
  # add to table
  overlap_ensembl_tb <- 
    tibble(coordinates = mcols(fpm_gr[queryHits(overlap_ensembl)])$coordinates,
           gene_id = names(exons_gr)[subjectHits(overlap_ensembl)]) %>% 
    dplyr::left_join(., genes_info, by = "gene_id") %>% 
    dplyr::select(coordinates, gene_hit = gene_id, gene_name, gene_biotype, gene_description) %>% 
    dplyr::mutate(gene_description = str_remove(gene_description, " \\[.*"))
  
  ### overlap with repeats
  overlap_rmsk <- findOverlaps(fpm_gr, rmsk_gr, ignore.strand = T)
  
  # add to table
  overlap_rmsk_tb <- tibble(coordinates = mcols(fpm_gr[queryHits(overlap_rmsk)])$coordinates,
                            repeat_hit_name = mcols(rmsk_gr[queryHits(overlap_rmsk)])$repName, 
                            repeat_hit_family = mcols(rmsk_gr[queryHits(overlap_rmsk)])$repFamily)
  
  
  ### add to table
  fpm_tb_annotated <- 
    fpm_tb %>% 
    dplyr::left_join(., overlap_ensembl_tb, by = "coordinates") %>% 
    dplyr::left_join(., overlap_rmsk_tb, by = "coordinates") %>% 
    dplyr::arrange(-FPM) 
  
  # save
  fpm_tb_annotated %T>% 
    readr::write_csv(., file.path(outpath, str_c("FPM.miR-205", "susScr11", name, "seed_match", "top100.annotated", "csv",  sep = ".")))
    
  # return 
  return(name)
  
})


