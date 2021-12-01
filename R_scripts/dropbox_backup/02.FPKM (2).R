### INFO: Count reads and get FPKM of genes
### DATE: Tue Dec 11 22:41:48 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression")

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
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.20181211.sampleTable.clean.csv"

# reduced exons path
exons_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.reducedExons.RDS"

# info about genes path
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

# summarizedExperiment path
se_path <- file.path(inpath, "Lnc1_KO.2018_Dec.GRCm38.91.reducedExons.summarizedOverlaps.RDS")

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path)

# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read summarizedExperiment
se <- readRDS(file = se_path) 

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.tibble(., rownames = "gene_id") %>%
  set_colnames(., str_remove(colnames(.), "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width),
                sample_id = str_remove(sample_id, "\\.SE")) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  tidyr::spread(key = sample_id, value = fpkm) %>%
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %T>%
  readr::write_csv(., path = file.path(outpath, "Lnc1_KO.2018_Dec.GRCm38.91.reducedExons.FPKM.csv"))
