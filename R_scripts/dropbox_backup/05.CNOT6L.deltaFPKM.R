### INFO: produce scatterplots of ratio of expression in Fugaku and CNOT6L data
### DATE: Sun Mar 11 00:36:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/fpkm")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(VennDiagram)
library(geneplotter)
library(GeneOverlap)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

# results path 
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# CNOT6L FPKM
fpkm_CNOT6L <- readr::read_csv(file = file.path(results_path, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# CNOT6L significantly diff. exp. genes
results_CNOT6L <- 
  list.files(file.path(results_path, "results"), "diffExp.CNOT6L.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# get genes which are significantly upregulated in MII
results_MII <- 
  results_CNOT6L %>% 
  dplyr::filter(comparison == "MII.KO_vs_WT", log2FoldChange > 0) %>% 
  dplyr::select(gene_id, gene_name, log2FoldChange, MII_WT, MII_KO) 

# save FPKM
results_MII %>% 
  dplyr::mutate(FPKM.delta = MII_KO - MII_WT) %>% 
  dplyr::rename(FPKM.MII_KO = MII_KO, FPKM.MII_WT = MII_WT) %T>%
  readr::write_csv(., path = file.path(outpath, "CNOT6L.signif_up.MII.KOvsWT.FPKM.csv"))

### calculate TPM
# total sum of FPKMs
fpkm_total <- 
  fpkm_CNOT6L %>% 
  dplyr::select(-gene_id) %>% 
  dplyr::summarise_all(.funs = funs(sum(.))) %>% 
  tidyr::gather(key = sample_id, value = fpkm_total) 

# divide FPKM by total FPKM in sample
tpm_CNOT6L <-
  fpkm_CNOT6L %>%
  tidyr::gather(key = sample_id, value = fpkm, -gene_id) %>%
  dplyr::left_join(., fpkm_total, by = "sample_id") %>%
  dplyr::mutate(tpm = (fpkm / fpkm_total) * 10^6) %>%
  dplyr::select(gene_id, sample_id, tpm) %>%
  tidyr::spread(key = sample_id, value = tpm) %T>%
  readr::write_csv(x = ., path = file.path(results_path, "ensembl.GRCm38.89.CNOT6L.avgTPM.csv"))

# get TPMs - divide FPKM by total FPKM
results_MII_tpm <-
  tpm_CNOT6L %>%
  dplyr::filter(gene_id %in% results_MII$gene_id) %>% 
  dplyr::select(gene_id, TPM.MII_KO = MII_KO, TPM.MII_WT = MII_WT) %>% 
  readr::write_csv(., path = file.path(outpath, "CNOT6L.signif_up.MII.KOvsWT.TPM.csv"))




