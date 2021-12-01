### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/tissues.Dicer_MT_HET.2021_Apr/Analysis/chr6.HET_vs_WT.expression")

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
genome_dir <- "/common/DB/genome_reference/rat/rn6.Rnor_6.0.GCA_000001895.4"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# expression path
expression_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/tissues.Dicer_MT_HET.2021_Apr/Analysis/expression"
expression_path <- list.files(expression_path, pattern = str_c("ensembl\\.", ensembl_version, ".*\\.FPKM\\.csv$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKM table
fpkm_tb <- readr::read_csv(expression_path)

######################################################## MAIN CODE
# filter FPKM values
fpkm_clean <- 
  fpkm_tb %>% 
  dplyr::filter(str_detect(coordinates, "^chr6")) %>% 
  dplyr::select(gene_id, starts_with("s_")) %>% 
  tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "fpkm") %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "MT_HET|WT"), 
                tissue = str_remove_all(sample_id, "^s_|_MT_HET.*|_WT.*")) %>% 
  dplyr::select(-sample_id) %>% 
  tidyr::pivot_wider(id_cols = c(gene_id, tissue), names_from = genotype, values_from = fpkm, names_prefix = "fpkm.") %>%
  dplyr::filter(fpkm.WT >= 1) %>% 
  dplyr::filter(fpkm.MT_HET <= (0.5*fpkm.WT)) %>% 
  tidyr::pivot_longer(-c(gene_id, tissue), names_to = "genotype", values_to = "fpkm", names_prefix = "fpkm.") %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "MT_HET"))) %>%
  dplyr::arrange(tissue, genotype) %>%  
  tidyr::pivot_wider(id_cols = gene_id, names_from = c(tissue, genotype), values_from = fpkm) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name, gene_description), by = "gene_id")
  
# save table
readr::write_csv(fpkm_clean, file.path(outpath, "ensembl.99.Rnor_6.0.chr6.HET_vs_WT.genes.csv"))




