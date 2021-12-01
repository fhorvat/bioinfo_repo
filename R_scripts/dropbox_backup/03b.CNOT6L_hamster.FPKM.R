### INFO: get expression of all genes in CNOT6L hamster data
### DATE: Sat Mar 10 20:38:05 2018
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
exons_path <- "/common/WORK/fhorvat/reference/hamster_golden/mesAur1/vfranke/ensembl.refSeq_annotation.MesAur1.0.91.20180309.reducedExons.RDS"

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1_vfranke"

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
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".total.bam", "")) %>% 
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id")

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "ensembl.mesAur1.CNOT6L.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.mesAur1.CNOT6L.summarizedOverlaps.RDS")) 

# ### FPKM
# # get data.frame of counts, transform to FPKM
# fpkm_df <-
#   assay(se) %>%
#   as.data.frame(.) %>%
#   tibble::rownames_to_column(., var = "gene_id") %>%
#   as.tibble(.) %>%
#   set_colnames(., str_replace(colnames(.), ".total.bam", "")) %>%
#   tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
#   dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
#   dplyr::left_join(., exons_width, by = "gene_id") %>%
#   dplyr::mutate(library_size = round(library_size / 1E6, 6),
#                 width = round(width / 1E3, 3),
#                 fpm = (counts / library_size),
#                 fpkm = (fpm / width)) %>%
#   dplyr::select(gene_id, sample_id, fpkm) %>%
#   dplyr::mutate(stage = str_replace_all(sample_id, "_r[1-3]{1}", "")) %>% 
#   dplyr::group_by(gene_id, stage) %>% 
#   dplyr::summarise(average_fpkm = mean(fpkm)) %>% 
#   dplyr::ungroup(.) %>% 
#   tidyr::spread(key = stage, value = average_fpkm) %T>%
#   readr::write_csv(., path = file.path(outpath, "ensembl.mesAur1.CNOT6L.avgFPKM.csv"))

### filter FPKM table
# read FPKM table
fpkm_df <- readr::read_csv(file = file.path(outpath, "ensembl.mesAur1.CNOT6L.avgFPKM.csv"))

# get FPKM of chosen genes
fpkm_filt <-
  fpkm_df %>%
  dplyr::right_join(., genes_info %>% dplyr::select(gene_id = mauratus_gene_id, gene_name), by = "gene_id") %>%
  dplyr::select(gene_name, everything(), -gene_id) %>%
  dplyr::arrange(gene_name) %T>%
  readr::write_csv(., path = file.path(outpath, "CNOT.mesAur1.CNOT6L.avgFPKM.csv"))

