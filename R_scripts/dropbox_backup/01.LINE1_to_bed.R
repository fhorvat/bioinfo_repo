### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/map_clusters_to_Siomi/LINE1_expression/LINE1_reads")

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

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# clean repeatMasker
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# read clean repeatMasker
rmsk_tb <- readr::read_delim(file = rmsk_clean_path, delim = "\t")

######################################################## MAIN CODE
### filter repeatMasker table
# get SINE, LINE and LTRs
rmsk_raw <- 
  rmsk_tb %>% 
  dplyr::filter(repClass %in% c("LINE"), 
                repFamily == "L1")

# format as table
rmsk_coords <- 
  rmsk_raw %>% 
  dplyr::mutate(gene_id = str_c(seqnames, start, end, rmsk_id, repName, sep = ".")) %>% 
  tidyr::unite(col = "coords", seqnames, start, end, sep = " ") %>% 
  dplyr::select(coords, strand, rmsk_id, gene_id, repName, repFamily, repClass)

# # save
# readr::write_csv(rmsk_coords, file.path(outpath, "rmsk.Siomi.20210519.LINE1.csv"))

### save as .bed
# to GRanges
rmsk_bed <- 
  rmsk_coords %>% 
  dplyr::select(coords, strand, rmsk_id, gene_id) %>% 
  tidyr::separate(coords, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# set names
names(rmsk_bed) <- mcols(rmsk_bed)$gene_id

# save
rtracklayer::export.bed(rmsk_bed, file.path(outpath, "rmsk.Siomi.20210519.LINE1.bed"))

