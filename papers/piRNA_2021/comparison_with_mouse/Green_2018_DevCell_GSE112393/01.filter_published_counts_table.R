### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Green_2018_DevCell_GSE112393/Analysis/expression")

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

library(Seurat)
library(patchwork)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw counts table
count_tb_path <- file.path(inpath, "GSE112393_MergedAdultMouseST25_DGE.txt.gz")

######################################################## READ DATA
# read raw published count table
count_tb <- 
  data.table::fread(count_tb_path) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "V1") %>% 
  as.matrix(.)

######################################################## MAIN CODE
# The starting pool of 59,313 cells and 37,241 genes was first selected by the cell size and integrity filters – 
# cells with ???500 detected genes per cell or with >= 10% of transcripts corresponding to mitochondria-encoded genes were removed. 
# We then removed low abundance genes that were detected in <= 15 cells or with <= 20 UMIs summed across all the retained cells. 
# These filters resulted in 
# 34,633 cells and 24,947 genes

### cell counts = columns
# for each cell count number of genes with at least one count
cell_genes_number <- colSums(count_tb > 0)

# for each cell calculate the percentage of reads coming from mitochondria-encoded genes
mito_reads <- colSums(count_tb[str_detect(rownames(count_tb), "^mt-"), ])
all_reads <- colSums(count_tb)
mito_percentage <- (mito_reads / all_reads) * 100

# get cells which have more than 500 detected genes AND which have less than 10% transcripts coming from mitochonria-encoded genes
cell_names <- colnames(count_tb)[cell_genes_number > 500 & mito_percentage < 10]

# filter table
count_tb.cells_filt <- count_tb[, cell_names]


### gene counts = rows
# for each gene count number of cells in which is detected
genes_cell_number <- rowSums(count_tb.cells_filt > 0)

# for each gene sum total number of reads
genes_count_sum <- rowSums(count_tb.cells_filt)

# get genes which are detected in more than 15 cells AND which have more than 20 total counts in all cells
gene_names <- rownames(count_tb.cells_filt)[genes_cell_number > 15 & genes_count_sum > 20] 

# filter genes
count_tb.cells_and_genes_filt <- count_tb.cells_filt[gene_names, ]

# This originally returns 34,633 cells and 24,923 genes, which is 24 genes fewer than they claim in the paper
# While checking code they published on GitHub, I realized they manually add 24 genes of interest which don't pass <= 20 UMIs / 15 cells per gene filtering. 
# These genes of interests are: 
gene_names_of_interest <- c("Ngf", "T", "Tbx20", "Tert", "Olig2", "Prl", "Hmx2", "Itgb6", "Mageb1", "Neurod1", "Neurod2",
                            "Ptchd4", "Prl6a1", "Fgf17", "Hpx", "Npy4r", "Fgf6", "Foxf2", "Prokr2", "Sry", "Ccl3", "Prf1", "Pax8", "Dcx",
                            "Ly6g5c", "Ereg", "Pgr", "Sycp2l", "Prnd", "Cysltr1", "Pdf")

# so I added this genes of interest as well
gene_names <- unique(c(gene_names, gene_names_of_interest))

# filter again, this is now returning 34,633 cells and 24,947 genes as stated in paper
count_tb.cells_and_genes_filt <- count_tb.cells_filt[gene_names, ]

# save filtered table
saveRDS(count_tb.cells_and_genes_filt, file.path(outpath, count_tb_path %>% basename(.) %>% str_replace(., "\\.txt\\.gz", ".filtered.matrix.RDS")))



