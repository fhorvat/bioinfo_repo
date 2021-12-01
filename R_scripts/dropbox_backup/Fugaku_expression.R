### INFO: 
### DATE: Mon Jun 17 17:11:29 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Analysis/expression")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# Fugaku summarizedOverlaps .RDS path
se_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.Fugaku.se.RDS"

# sample table path 
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Data/Documentation/Fugaku.20190426.sampleTable.csv"

######################################################## READ DATA
# read genes info
genes_info <- data.table::fread(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read summarizedOverlaps
se <- readRDS(se_path)

# read sample table
sample_table <- data.table::fread(sample_table_path)

######################################################## MAIN CODE
# get total length of all exons for each gene
exons_table <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .) %>% 
  as.data.table(.)

# FPKM
fpkm_long_tb <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[, sample_id := str_remove(sample_id, ".genome.Aligned.sortedByCoord.out.bam|.total.bam")] %>% 
  .[sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)), 
                                                        stage = i.stage)] %>% 
  .[exons_table, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3))] %>% 
  .[]

# wide and save
fpkm_wide_tb <- 
  fpkm_long_tb %>% 
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %>% 
  .[genes_info, on = "gene_id", `:=`(gene_symbol = i.gene_name, 
                                     coordinates_mm10 = str_c(i.seqnames, i.start, i.end, sep = " "))] %>% 
  as_tibble(.) %>%  
  dplyr::select(ensembl_id = gene_id, gene_symbol, coordinates_mm10, s_GV.WE.PE, s_MII.WE.PE, s_1cell.WE.PE, s_2cell.WE.PE, s_4cell.WE.PE, s_Molura.WE.PE, s_Blast.WE.PE) %>% 
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".wide.FPKM.csv")))
