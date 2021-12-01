### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Pal_2020_unpub_GSE140545/Analysis/expression")

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

library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/human/hg38.GRCh38.GCA_000001405.15"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.reducedExons\\.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(inpath, ".*\\.FPKM\\.csv")

# bigwig track path
mapped_path <- file.path(inpath, "../../Data/Mapped/STAR_hg38")
bw_path <- list.files(mapped_path, "s_RPE_r1\\.PE.*\\.bw", full.names = T)

####################################################### READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

# read bigWig track
coverage_gr <- rtracklayer::import(bw_path)

######################################################## MAIN CODE
# overlap coverage and genes
overlaps <- findOverlaps(exons_gr, coverage_gr, ignore.strand = T)

# extract overlaps
gene_ids_list <- names(exons_gr[queryHits(overlaps)])
coverage_list <- mcols(coverage_gr[subjectHits(overlaps)])$score



# overlap coverage and genes
overlaps_test <- findOverlaps(exons_gr$ENSG00000182492, coverage_gr, ignore.strand = T)

# extract overlaps
coverage_list_test <- coverage_gr[subjectHits(overlaps_test)]


# get table
gene_coverage <- 
  tibble(gene_id = gene_ids_list, 
         rpm = coverage_list) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::slice_max(order_by = rpm, n = 1, with_ties = F) %>% 
  dplyr::ungroup(.) 

# get log2FC between RPE and RPE progenitors
fpkm_logfc <- 
  fpkm_tb %>% 
  dplyr::select(-s_photoreceptor_precursors_r1.PE) %>% 
  dplyr::mutate(logFC.RPE_vs_RPEpro = log2(s_RPE_r1.PE / s_RPE_progenitors_r1.PE)) %>% 
  dplyr::left_join(., gene_coverage, by = "gene_id") %>% 
  dplyr::select(gene_id, gene_name, s_RPE_progenitors_r1.PE:s_retinal_progenitors_r1.PE, 
                logFC.RPE_vs_RPEpro, max_CPM_peak_in_RPEs = rpm, coordinates:gene_description) %>% 
  dplyr::mutate(max_CPM_peak_in_RPEs = replace(max_CPM_peak_in_RPEs, is.na(max_CPM_peak_in_RPEs), 0)) 

# save table
readr::write_csv(fpkm_logfc, file.path(outpath, "ensembl.99.FPKM.logFC_RPEs_vs_RPEprogenitors.csv"))

# filter
fpkm_filt <- 
  fpkm_logfc %>% 
  dplyr::filter(logFC.RPE_vs_RPEpro > 5) %>% 
  dplyr::filter(max_CPM_peak_in_RPEs > 25) %>% 
  dplyr::arrange(-logFC.RPE_vs_RPEpro)

# save
readr::write_csv(fpkm_filt, file.path(outpath, "ensembl.99.FPKM.logFC_RPEs_vs_RPEprogenitors.filtered.csv"))



