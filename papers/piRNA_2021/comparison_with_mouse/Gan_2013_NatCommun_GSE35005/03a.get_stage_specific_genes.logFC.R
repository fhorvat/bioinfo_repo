### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Gan_2013_NatCommun_GSE35005")

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

library(openxlsx)

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
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# mouse expression path
mouse_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Gan_2013_NatCommun_GSE35005/Analysis/expression"
mouse_fpkm_path <- list.files(mouse_path, str_c("ensembl.", ensembl_version, ".*\\FPKM_mean.csv"), full.names = T)
mouse_results_path <- list.files(file.path(mouse_path, "results.ensembl.99.GRCm38.p6.20200415.UCSCseqnames"), 
                                 str_c("protein_coding\\.diffExp\\.DESeq2.*\\.all_results\\.xlsx"), full.names = T)

# hamster-mouse orthologs path
ortho_path <- file.path(inpath, "..", "ensembl.99.mouse_vs_goldHamster.one2one_homologs.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKM table
fpkm_tb <- readr::read_csv(mouse_fpkm_path)

# read orthologs
ortho_tb <- readr::read_csv(ortho_path)

### read mouse DESeq2 results
# get sheet names
sheet_names <- openxlsx::getSheetNames(mouse_results_path)

# loop through sheet names
mouse_results_list <- purrr::map(sheet_names, function(sheet_name){
  
  # read one table, convert to tibble, select columns
  mouse_tb <-
    openxlsx::read.xlsx(mouse_results_path, sheet = sheet_name) %>% 
    as_tibble(.) %>% 
    dplyr::select(gene_id, log2FoldChange) %>% 
    magrittr::set_colnames(., c("gene_id", str_c("logFC.", sheet_name)))

  # return
  return(mouse_tb)
  
})

######################################################## MAIN CODE
# join results into one table
mouse_tb <- 
  purrr::reduce(mouse_results_list, left_join, by = "gene_id") %>% 
  dplyr::left_join(., fpkm_tb %>% dplyr::select(gene_id, SG_A:round_ST), by = "gene_id")

### filter stage specific genes
# set FPKM and log2FC cutoff
fpkm_cutoff <- 5
log2fc_cutoff <- 2
  
# stages = primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC, round_ST, elongative_ST
# primitive spermatogonia A (primitive_SG_A)
primitive_SG_A <- 
  mouse_tb %>% 
  dplyr::filter(primitive_SG_A > fpkm_cutoff) %>% 
  dplyr::filter(logFC.SG_A_vs_primitive_SG_A < 0) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "primitive_SG_A")

# spermatogonia A (SG_A)
SG_A <- 
  mouse_tb %>% 
  dplyr::filter(SG_A > fpkm_cutoff) %>% 
  dplyr::filter(logFC.SG_A_vs_primitive_SG_A > log2fc_cutoff) %>% 
  dplyr::filter(logFC.SG_B_vs_SG_A <= 0)  %>%
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "SG_A")

# spermatogonia B (SG_B)
SG_B <- 
  mouse_tb %>% 
  dplyr::filter(SG_B > fpkm_cutoff) %>% 
  dplyr::filter(logFC.SG_B_vs_SG_A > log2fc_cutoff) %>% 
  dplyr::filter(logFC.leptotene_SC_vs_SG_B <= 0) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "SG_B")

# leptotene spermatocyes (leptotene_SC)
leptotene_SC <- 
  mouse_tb %>% 
  dplyr::filter(leptotene_SC > fpkm_cutoff) %>% 
  dplyr::filter(logFC.leptotene_SC_vs_SG_B > log2fc_cutoff) %>% 
  dplyr::filter(logFC.pachytene_SC_vs_leptotene_SC <= 0) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "leptotene_SC")

# pachytene spermatocyes (pachytene_SC)
pachytene_SC <- 
  mouse_tb %>% 
  dplyr::filter(pachytene_SC > fpkm_cutoff) %>% 
  dplyr::filter(logFC.pachytene_SC_vs_leptotene_SC > log2fc_cutoff) %>% 
  dplyr::filter(logFC.round_ST_vs_pachytene_SC <= 0) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "pachytene_SC")

# round spermatids (round_ST)
round_ST <- 
  mouse_tb %>% 
  dplyr::filter(round_ST > fpkm_cutoff) %>% 
  dplyr::filter(logFC.round_ST_vs_pachytene_SC > log2fc_cutoff) %>% 
  dplyr::filter(logFC.elongative_ST_vs_round_ST <= 0) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "round_ST") 
  
# elongative spermatids (elongative_ST)
elongative_ST <- 
  mouse_tb %>% 
  dplyr::filter(elongative_ST > fpkm_cutoff) %>%
  dplyr::filter(logFC.elongative_ST_vs_round_ST > log2fc_cutoff) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::mutate(stage = "elongative_ST")


### join all, find hamster orthologs
stage_specific_genes <- 
  list(primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC, round_ST, elongative_ST) %>% 
  bind_rows(.) %>% 
  dplyr::filter(!duplicated(gene_id)) %>% 
  dplyr::mutate(stage = factor(stage, levels = c("primitive_SG_A", "SG_A", "SG_B", "leptotene_SC", "pachytene_SC", "round_ST", "elongative_ST"))) %>% 
  left_join(., ortho_tb %>% dplyr::select(gene_id = mouse_gene_id, hamster_gene_id), by = "gene_id") %>% 
  dplyr::filter(!is.na(hamster_gene_id)) %T>%
  readr::write_csv(., file.path(outpath, "Gan_2013_NatCommun_GSE35005.stage_specific_genes.hamster_gene_id.csv"))


  