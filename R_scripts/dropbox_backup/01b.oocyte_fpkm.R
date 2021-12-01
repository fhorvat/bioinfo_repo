### INFO: expression of lnc1 and other genes in ENCODE dataset
### DATE: Tue Aug 21 14:25:50 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/lncRNA_expression/lnc1_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path 
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
gene_info_path <- list.files(genome_path, "ensembl.91.GRCm38.p5.*.UCSCseqnames.geneInfo.csv", full.names = T)

# path to ENCODE mapped samples
samples_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417"
  
# sample table path
sample_table_path <- list.files(file.path(samples_path, "Analysis"), pattern = ".*.sample_table.csv", full.names = T)

# Fugaku FPKM expression table
fugaku_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/expression/ensembl.91.GRCm38.p5.20180512.rmskFiltered.Fugaku.fpkm.csv"

# CNOT6L FPKM expression table
cnot6l_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read Fugaku FPKM table
fugaku_fpkm <- readr::read_csv(fugaku_path)

######################################################## MAIN CODE
# get topN expressed genes in Fugaku's GV
genes_top <- 
  fugaku_fpkm %>% 
  dplyr::select(gene_id, s_GV.WE.PE) %>% 
  dplyr::left_join(., gene_info, by = "gene_id") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lincRNA")) %>%
  dplyr::arrange(desc(s_GV.WE.PE)) %>% 
  dplyr::slice(1:20) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(gene_id, gene_name, fpkm = s_GV.WE.PE, coordinates, gene_biotype, gene_description) %>% 
  dplyr::mutate(gene_name = factor(gene_name, levels = gene_name)) %T>%
  readr::write_csv(., path = file.path(outpath, "fugaku_GV.WE.top20.fpkm.csv"))

# visualize
ggplot(genes_top, aes(x = gene_name, y = fpkm, fill = gene_biotype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red3", "gray30")) +
  ylab("Fugaku GV FPKM") + 
  xlab("gene") + 
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, "fugaku_GV.WE.barplot.top20.png"), width = 12, height = 10)


