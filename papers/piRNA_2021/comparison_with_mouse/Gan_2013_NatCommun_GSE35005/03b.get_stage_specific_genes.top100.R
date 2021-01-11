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

# hamster expression path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression"
hamster_tb_path <- list.files(hamster_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKM table
fpkm_tb <- readr::read_csv(mouse_fpkm_path)

# read orthologs
ortho_tb <- readr::read_csv(ortho_path)

# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

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

# stages = primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC, round_ST, elongative_ST
# primitive spermatogonia A (primitive_SG_A)
primitive_SG_A_vs_SG_A <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(primitive_SG_A, SG_A), any_vars(. > fpkm_cutoff))
primitive_SG_A_vs_SG_A <-  
  bind_rows(dplyr::top_n(primitive_SG_A_vs_SG_A, 50, logFC.SG_A_vs_primitive_SG_A), 
            dplyr::top_n(primitive_SG_A_vs_SG_A, 50, -logFC.SG_A_vs_primitive_SG_A)) %>% 
  dplyr::select(gene_id, logFC = logFC.SG_A_vs_primitive_SG_A) %>% 
  dplyr::mutate(stage = "primitive_SG_A_vs_SG_A")

# spermatogonia A (SG_A)
SG_A_vs_SG_B <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(SG_A, SG_B), any_vars(. > fpkm_cutoff))
SG_A_vs_SG_B <-  
  bind_rows(dplyr::top_n(SG_A_vs_SG_B, 50, logFC.SG_B_vs_SG_A), 
            dplyr::top_n(SG_A_vs_SG_B, 50, -logFC.SG_B_vs_SG_A)) %>% 
  dplyr::select(gene_id, logFC = logFC.SG_B_vs_SG_A) %>% 
  dplyr::mutate(stage = "SG_A_vs_SG_B")

# spermatogonia B (SG_B)
SG_B_vs_leptotene_SC <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(SG_B, leptotene_SC), any_vars(. > fpkm_cutoff))
SG_B_vs_leptotene_SC <-  
  bind_rows(dplyr::top_n(SG_B_vs_leptotene_SC, 50, logFC.leptotene_SC_vs_SG_B), 
            dplyr::top_n(SG_B_vs_leptotene_SC, 50, -logFC.leptotene_SC_vs_SG_B)) %>% 
  dplyr::select(gene_id, logFC = logFC.leptotene_SC_vs_SG_B) %>% 
  dplyr::mutate(stage = "SG_B_vs_leptotene_SC")

# leptotene spermatocyes (leptotene_SC)
leptotene_SC_vs_pachytene_SC <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(leptotene_SC, pachytene_SC), any_vars(. > fpkm_cutoff))
leptotene_SC_vs_pachytene_SC <-  
  bind_rows(dplyr::top_n(leptotene_SC_vs_pachytene_SC, 50, logFC.pachytene_SC_vs_leptotene_SC), 
            dplyr::top_n(leptotene_SC_vs_pachytene_SC, 50, -logFC.pachytene_SC_vs_leptotene_SC)) %>% 
  dplyr::select(gene_id, logFC = logFC.pachytene_SC_vs_leptotene_SC) %>% 
  dplyr::mutate(stage = "leptotene_SC_vs_pachytene_SC")

# pachytene spermatocyes (pachytene_SC)
pachytene_SC_vs_round_ST <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(pachytene_SC, round_ST), any_vars(. > fpkm_cutoff))
pachytene_SC_vs_round_ST <-  
  bind_rows(dplyr::top_n(pachytene_SC_vs_round_ST, 50, logFC.round_ST_vs_pachytene_SC), 
            dplyr::top_n(pachytene_SC_vs_round_ST, 50, -logFC.round_ST_vs_pachytene_SC)) %>% 
  dplyr::select(gene_id, logFC = logFC.round_ST_vs_pachytene_SC) %>% 
  dplyr::mutate(stage = "pachytene_SC_vs_round_ST")

# round spermatids (round_ST)
round_ST_vs_elongative_ST <- 
  mouse_tb %>% 
  dplyr::filter_at(.vars = vars(round_ST, elongative_ST), any_vars(. > fpkm_cutoff))
round_ST_vs_elongative_ST <-  
  bind_rows(dplyr::top_n(round_ST_vs_elongative_ST, 50, logFC.elongative_ST_vs_round_ST), 
            dplyr::top_n(round_ST_vs_elongative_ST, 50, -logFC.elongative_ST_vs_round_ST)) %>% 
  dplyr::select(gene_id, logFC = logFC.elongative_ST_vs_round_ST) %>% 
  dplyr::mutate(stage = "round_ST_vs_elongative_ST")


### join all, find hamster orthologs
stage_specific_genes <- 
  list(primitive_SG_A_vs_SG_A, SG_A_vs_SG_B, SG_B_vs_leptotene_SC, leptotene_SC_vs_pachytene_SC, pachytene_SC_vs_round_ST, round_ST_vs_elongative_ST) %>% 
  bind_rows(.) %>% 
  # dplyr::filter(!duplicated(gene_id)) %>%
  dplyr::mutate(stage = factor(stage, levels = c("primitive_SG_A_vs_SG_A", "SG_A_vs_SG_B", "SG_B_vs_leptotene_SC", "leptotene_SC_vs_pachytene_SC", 
                                                 "pachytene_SC_vs_round_ST", "round_ST_vs_elongative_ST"))) %>% 
  left_join(., ortho_tb %>% dplyr::select(gene_id = mouse_gene_id, hamster_gene_id), by = "gene_id") %>% 
  dplyr::filter(!is.na(hamster_gene_id)) %>%
  dplyr::select(hamster_gene_id, mouse_gene_id = gene_id, comparison = stage, mouse_logFC = logFC) %T>% 
  readr::write_csv(., file.path(outpath, "Gan_2013_NatCommun_GSE35005.stage_specific_genes.hamster_gene_id.top_100.csv"))

# join with FPKM values, save
stage_specific_genes %>% 
  dplyr::left_join(., fpkm_tb %>% dplyr::select(gene_id, primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC, round_ST, elongative_ST),
                   by = c("mouse_gene_id" = "gene_id")) %>% 
  dplyr::left_join(., hamster_tb, by = c("hamster_gene_id" = "gene_id")) %>% 
  dplyr::arrange(comparison, desc(mouse_logFC)) %T>%
  readr::write_csv(., file.path(outpath, "Gan_2013_NatCommun_GSE35005.stage_specific_genes.hamster_gene_id.top_100.FPKM.csv"))


