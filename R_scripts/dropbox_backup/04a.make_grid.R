### INFO: 
### DATE: Thu Oct 25 10:46:52 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.geneInfo.csv$", full.names = T)

# ENCODE FPKM path
encode_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.FPKM.csv"

# Fugaku FPKM path
fugaku_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.Fugaku.FPKM.csv"

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read info about chosen genes
genes_info <- data.table::fread(file = genes_info_path)

# read ENCODE FPKM table
encode_fpkm <- data.table::fread(file = encode_path)

# read Fugaku FPKM table
fugaku_fpkm <- data.table::fread(file = fugaku_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

##### ENCODE ####
# get FPKM's for only protein coding genes, get mean expression across all tissue
encode_fpkm_clean <- 
  encode_fpkm[(gene_id %in% protein_genes) & (sample_id != "ovary")] %>% 
  .[, list(fpkm_tissue = round(mean(fpkm), 3)), by = gene_id]

### Fugaku ####
# get FPKM's for only protein coding genes in oocyte, bin FPKM to 100 bins of equal size
fugaku_fpkm_clean <- 
  fugaku_fpkm[(gene_id %in% protein_genes) & (sample_id == "s_GV.WE") & (fpkm > 0), -"sample_id"] %>% 
  .[encode_fpkm_clean, on = "gene_id", fpkm_tissue := i.fpkm_tissue] %>% 
  .[, fpkm_relative := (fpkm - fpkm_tissue)] %>% 
  .[order(fpkm), bin_fpkm := as.numeric(Hmisc::cut2(fpkm, g = 100))] %>% 
  .[order(fpkm_relative), bin_fpkm_relative := as.numeric(Hmisc::cut2(fpkm_relative, g = 100))] %>% 
  .[]

# save binned FPKMs
readr::write_csv(fugaku_fpkm_clean, file.path(outpath, "grid.Fugaku.Encode.binned_fpkm.100.csv"))
