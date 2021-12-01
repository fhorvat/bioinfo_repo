### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/compare_mapping")

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
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# expression path
STAR_path <- file.path(inpath, "../expression")
bbmap_path <- file.path(inpath, "../expression.bbmap")

# FPKM path
STAR_fpkm_path <- list.files(STAR_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_long\\.csv"), full.names = T)
bbmap_fpkm_path <- list.files(bbmap_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_long\\.csv"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKM tables
STAR_fpkm <- readr::read_csv(STAR_fpkm_path)
bbmap_fpkm <- readr::read_csv(bbmap_fpkm_path)

######################################################## MAIN CODE
# join two tables, plot scatter plot
fpkm_all <- 
  STAR_fpkm %>% 
  dplyr::select(sample_id, gene_id, STAR = fpkm) %>% 
  dplyr::left_join(., bbmap_fpkm %>% dplyr::select(sample_id, gene_id, bbmap = fpkm), by = c("sample_id", "gene_id")) %>% 
  
  dplyr::filter(2*STAR < bbmap) %>% 
  dplyr::left_join(., genes_info, by = "gene_id")
  
  dplyr::mutate_at(vars(STAR, bbmap), ~(log2(. + 0.1)))

# plot 
scatter_plot <- 
  ggplot(data = fpkm_all, aes(x = STAR, y = bbmap)) + 
  geom_abline(intercept = 0, color = "red", size = 1.5) +
  geom_point(size = 1) +
  facet_wrap(vars(sample_id), nrow = 4, ncol = 3) + 
  xlab("STAR") + 
  ylab("bbmap") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave(filename = file.path(outpath, "scatter_plot.compare_STAR_bbmap.FPKM.png"), width = 12, height = 12)

