### INFO: 
### DATE: Thu Oct 25 10:46:52 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/NGS_data")

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
encode_fpkm_clean <- encode_fpkm[(gene_id %in% protein_genes) & (sample_id != "ovary")]

# get FPKM's for only protein coding genes in oocyte, bin FPKM to 61 bins of equal size
fugaku_fpkm_relative <- 
  fugaku_fpkm[(gene_id %in% protein_genes) & (sample_id == "s_GV.WE")] %>% 
  .[, sample_id := "oocyte"] %>% 
  rbind(., encode_fpkm_clean) %>% 
  .[order(gene_id, fpkm)] %>% 
  .[, position := order(fpkm), by = "gene_id"] %>% 
  .[sample_id == "oocyte", ] %>%
  .[order(position), bin_relative := as.numeric(Hmisc::cut2(1:.N, g = 21))] %>% 
  .[order(position)]

### Fugaku ####
# get FPKM's for only protein coding genes in oocyte, bin FPKM to 61 bins of equal size
fugaku_fpkm_absolute <- 
  fugaku_fpkm[(gene_id %in% protein_genes) & (sample_id == "s_GV.WE") & (fpkm > 0), -"sample_id"] %>% 
  .[order(fpkm), bin_absolute := as.numeric(Hmisc::cut2(fpkm, g = 21))] %>% 
  .[fugaku_fpkm_relative, on = "gene_id", `:=`(bin_relative = i.bin_relative)] %>% 
  .[]

# save binned FPKMs
readr::write_csv(fugaku_fpkm_absolute, file.path(outpath, "grid.Fugaku.Encode.oocyte.21.bins.csv"))
