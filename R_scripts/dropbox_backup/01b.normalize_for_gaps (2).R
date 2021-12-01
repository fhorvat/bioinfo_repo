### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/count_N_in_genome")

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
library(Biostrings)
library(BSgenome.Maur.UCSC.MesAur1)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# RPM table path
rpm_tb_path <- file.path(inpath, "../recalculate_cluster_expression")
rpm_tb_path <- list.files(rpm_tb_path, ".*\\.FPM_mean\\.csv$", full.names = T, recursive = T)

# fasta path
fasta_path <- "../small_RNAseq"
fasta_path <- list.files(fasta_path, "\\.fa$", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE# set name
# separate for PIWIL1 and PIWIL3 clusters
purrr::map(c("PIWIL1", "PIWIL3"), function(piwil){
  
  # read RPM table
  rpm_tb <- readr::read_csv(rpm_tb_path[str_detect(rpm_tb_path, piwil)])
  
  # read fasta
  tile_seq <- Biostrings::readDNAStringSet(fasta_path[str_detect(fasta_path, piwil)])
  
  # set table name
  table_name <- 
    rpm_tb_path[str_detect(rpm_tb_path, piwil)] %>% 
    basename(.) 
  
  # count N's in fasta
  count_n <- vcountPattern(pattern = "N", tile_seq)
  
  # create table 
  nucleotide_tb <- tibble(gene_id = names(tile_seq), 
                          count_N = count_n)
  
  # get coordinates from RPM table
  tile_tb <- 
    rpm_tb %>% 
    dplyr::select(coordinates, gene_id) %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
    GRanges(.) %>% 
    as_tibble(.) %>% 
    dplyr::left_join(., nucleotide_tb, by = "gene_id") %>% 
    dplyr::mutate(width_norm = width - count_N) %>% 
    dplyr::select(gene_id, width = width_norm) 
  
  # save
  readr::write_csv(tile_tb, file.path(inpath, table_name %>% str_replace(., "FPM_mean", "norm_width")))
  
  # # normalize rpm
  # rpkm_tb <-
  #   rpm_tb %>%
  #   tidyr::pivot_longer(-c(gene_id, coordinates), values_to = "rpm", names_to = "sample_id") %>%
  #   dplyr::left_join(., tile_tb, by = "gene_id") %>%
  #   dplyr::filter(width > 0) %>%
  #   dplyr::mutate(width = (width / 1000),
  #                 rpkm = round((rpm / width), 3)) %>%
  #   tidyr::pivot_wider(id_cols = c(gene_id, coordinates), names_from = "sample_id", values_from = "rpkm")
  # 
  # # save
  # readr::write_csv(rpkm_tb, file.path(inpath, table_name %>% str_replace(., "FPM", "RPKM")))
  
  # return
  return(piwil)
  
})

