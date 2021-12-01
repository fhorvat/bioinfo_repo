### INFO: 
### DATE: Thu Oct 25 10:46:52 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/last_exon.let7_counts")

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

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.geneInfo.csv$", full.names = T)

# CNOT6L FPKM path
expression_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/datasets/CNOT6L/expression/ensembl.93.mouse_B6.mm10.FPKM.csv"

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Documentation/CNOT6L.sample_table.csv"

######################################################## READ DATA
# read info about chosen genes
genes_info <- data.table::fread(file = genes_info_path)

# read CNOT6L FPKM table
expression_fpkm <- data.table::fread(file = expression_path)

# read sample table
sample_table <- data.table::fread(file = sample_table_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

##### CNOT6L ####
# get FPKM's for only protein coding genes, get mean expression across all tissue
expression_fpkm_clean <- expression_fpkm[(gene_id %in% protein_genes)]

# get FPKM's for only protein coding genes in oocyte, get average FPKM, filter
expression_fpkm_avg <- 
  expression_fpkm_clean[(gene_id %in% protein_genes)] %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "fpkm") %>% 
  .[sample_table[, -"bam_path"], on = "sample_id", `:=`(genotype = i.genotype, 
                                                        stage = i.stage)] %>% 
  .[, list(avg_fpkm = round(mean(fpkm), 3)), by = .(gene_id, stage, genotype)] %>% 
  dcast(., gene_id ~ stage + genotype, value.var = "avg_fpkm") %>% 
  .[GV_WT >= 0.1, list(gene_id, GV_WT)]

# save
saveRDS(expression_fpkm_avg, file.path(outpath, "ensembl.93.mouse_B6.mm10.avgFPKM_above_0.1.RDS"))
