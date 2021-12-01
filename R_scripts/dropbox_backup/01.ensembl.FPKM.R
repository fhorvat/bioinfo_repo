#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get expression in mouse dataset
### DATE: Sun Jun 24 16:14:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/expression")

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

# set experiment name
experiment <- "Fugaku"
experiment_name <- "Fugaku"

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.91.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# # filtered reduced exons path
# exons_filt_path <- list.files(path = genome_dir, pattern = "ensembl.91.*UCSCseqnames.reducedExons.rmskFiltered.RDS$", full.names = T)

# filtered summarizedOverlaps path
se_filt_path <- list.files(inpath, "*.rmskFiltered.Fugaku.se.RDS", full.names = T)

# stats path
stats_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_new/log.Fugaku.stats_and_tracks.csv"
  
######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# # read ENSEMBL filtered reduced exons
# exons_filt_gr <- readRDS(file = exons_filt_path)

# read summarizedExperiment from RDS file
se_filt <- readRDS(file = se_filt_path)

# read stats and tracks
stats_df <- readr::read_csv(stats_path)

######################################################## MAIN CODE
# create sample table
sample_table <- 
  stats_df %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA)
  
# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

### FPKM
# get data.frame of counts, transform to FPKM, write
fpkm_df <-
  assay(se_filt) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  tibble::as.tibble(.) %>%
  magrittr::set_colnames(., stringr::str_remove(colnames(.), ".Aligned.sortedByCoord.out.bam|.total.bam")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  tidyr::spread(key = sample_id, value = fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, basename(se_filt_path) %>% 
                                         str_remove(., ".UCSCseqnames") %>% 
                                         str_replace(".se.RDS", ".fpkm.csv")))
