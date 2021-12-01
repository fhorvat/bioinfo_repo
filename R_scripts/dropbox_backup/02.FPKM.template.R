#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: count reads over ENSEMBL
### DATE: Wed Sep 26 17:40:45 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/ensembl_counts/%EXPERIMENT")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genomes/%EXPERIMENT"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.93.*reducedExons.RDS$", full.names = T)

# read additional info about genes
genes_info <- list.files(path = genome_dir, pattern = "ensembl.93.*geneInfo.csv$", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/%EXPERIMENT"

# get list of bam files
bam_path <- list.files(path = mapped_path, pattern = "*.total.bam$", full.names = T)

# stats&tracks path
stats_path <- list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T)

# summarizedExpriment path
se_path <- list.files(path = inpath, pattern = "ensembl.93.*se.RDS$", full.names = T)

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read stats and tracks table 
stats_df <- readr::read_csv(file = stats_path)

# read summarizedExpriment
se <- readRDS(file = se_path)

######################################################## MAIN CODE
# get library size
library_size_df <- 
  stats_df %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA)

# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_remove(colnames(.), ".total.bam")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., library_size_df, by = "sample_id") %>% 
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, ".PE")) %>%
  tidyr::spread(key = sample_id, value = fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("ensembl.93.%EXPERIMENT.FPKM.csv")))
