### INFO: 
### DATE: Fri Dec 11 17:05:38 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq/Analysis/expression.added_PIWIL3.stranded")

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

# differently expressed genes in testis path
diffExp_testis_path <- file.path(inpath, "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq", 
                                 "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx")

# differently expressed genes in oocyte path
diffExp_oocyte_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression.added_PIWIL3"
diffExp_oocyte_path <- file.path(diffExp_oocyte_path, "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq", 
                                 "protein_coding.diffExp.DESeq2.genotype.significant_results.xlsx")

######################################################## READ DATA
# read differently expressed genes
diffExp_testis <- openxlsx::read.xlsx(diffExp_testis_path) %>% as_tibble(.)
diffExp_oocyte <- openxlsx::read.xlsx(diffExp_oocyte_path) %>% as_tibble(.)

######################################################## MAIN CODE
# clean testis table
diffExp_testis_tb <- 
  diffExp_testis %>% 
  dplyr::select(gene_id, baseMean.testis = baseMean, log2FoldChange.testis = log2FoldChange, 
                Mov10l1_KO_8.5dpp.testis = Mov10l1_KO_8.5dpp.FPKM, Mov10l1_WT_8.5dpp.testis = Mov10l1_WT_8.5dpp.FPKM, 
                gene_name:gene_description)

# clean oocyte table
diffExp_oocyte_tb <- 
  diffExp_oocyte %>% 
  dplyr::select(gene_id, baseMean.GV = baseMean, log2FoldChange.GV = log2FoldChange, 
                Mov10l1_KO.GV = Mov10l_KO.FPKM, Mov10l1_WT.GV = Mov10l_WT.FPKM)

# join together to see what's common
diffExp_common <- 
  dplyr::inner_join(diffExp_testis_tb, diffExp_oocyte_tb, by = "gene_id") %>% 
  dplyr::select(gene_id, gene_name, 
                log2FoldChange.testis, log2FoldChange.GV, 
                baseMean_DESeq.testis = baseMean.testis, baseMean_DESeq.GV = baseMean.GV, 
                Mov10l1_WT_8.5dpp.testis, Mov10l1_KO_8.5dpp.testis, 
                Mov10l1_WT.GV, Mov10l1_KO.GV, 
                coordinates, strand,
                gene_biotype, gene_description) %T>% 
  readr::write_csv(., file.path(outpath, "protein_coding.diffExp.DESeq2..significant_results.shared.testis_8.5dpp.oocyte.xlsx"))
