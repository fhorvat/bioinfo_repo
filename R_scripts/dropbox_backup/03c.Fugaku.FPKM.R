### INFO: get expression of all genes in Fugaku's data
### DATE: Fri Mar 09 22:34:51 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

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

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reduced exons path
exons_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.reducedExons.RDS"

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_noMultimapFilter"

# info about genes path
genes_info_path <- file.path(outpath, "CNOT.GRCm38.89.20180305.geneInfo.csv")

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read stats and tracks table
stats_df <- readr::read_csv(file = list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T))

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# construct sample table
sample_table <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".Aligned.sortedByCoord.out.bam", "")) %>% 
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, mapped_genome_except_rDNA), by = "sample_id") %>% 
  dplyr::rename(library_size = mapped_genome_except_rDNA) %T>% 
  readr::write_csv(., path = file.path(outpath, "Fugaku.sample_table.csv"))

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "ensembl.GRCm38.89.Fugaku.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.GRCm38.89.Fugaku.summarizedOverlaps.RDS"))

# read sample table
sample_table <- readr::read_csv(file = file.path(outpath, "Fugaku.sample_table.csv"))

### FPKM
# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_", "") %>% 
                  str_replace(., "cell", "C")) %>%
  tidyr::spread(key = sample_id, value = fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, "ensembl.GRCm38.89.Fugaku.FPKM.csv"))

# get FPKM of chosen genes
fpkm_df %>%
  dplyr::right_join(., genes_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::select(gene_name, everything(), -gene_id) %>%
  dplyr::arrange(gene_name) %T>%
  readr::write_csv(., path = file.path(outpath, "CNOT.GRCm38.89.Fugaku.FPKM.csv"))
