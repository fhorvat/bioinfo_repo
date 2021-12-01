### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/review/GO_terms")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(goseq)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "hamster_oocyte_Mov10l.RNAseq"
experiment_name <- "hamster_oocyte_Mov10l.RNAseq"


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### results path
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", experiment)

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression.added_PIWIL3")
analysis_path <- file.path(analysis_path, "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq")

# differential expression results path
diffExp_path <- list.files(analysis_path, "protein_coding.diffExp.DESeq2.genotype.significant_results.xlsx", full.names = T)
diffExp_path_all <- list.files(analysis_path, "protein_coding.diffExp.DESeq2.genotype.all_results.xlsx", full.names = T)


### genome path
# set ensembl version
ensembl_version <- 99

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# genes info path
genes_info_path <- list.files(path = genome_path, 
                              pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), 
                              full.names = T)

# mouse homologs path
mouse_homologs_path <- list.files(genome_path, pattern = str_c("ensembl.", ensembl_version, ".*mmusculus_homologs.csv$"), full.names = T)

# mouse genes info
ensembl_version <- 99
genes_info_path_mouse <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
genes_info_path_mouse <- list.files(path = genes_info_path_mouse, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read gene info
genes_info <- readr::read_csv(genes_info_path)

# read mouse homologs
mouse_homologs <- readr::read_csv(mouse_homologs_path)

# read mouse genes info
genes_info_mouse <- readr::read_csv(genes_info_path_mouse)


### read list of diff. expressed genes
# set comparison
comparison_list <- c("Mov10l_KO_vs_Mov10l_WT", "Mov10l_HET_vs_Mov10l_WT")

# read different comparisons - significant genes
diffExp_results <- purrr::map(comparison_list, function(comparison){
  
  # read sheet
  read.xlsx(xlsxFile = diffExp_path, sheet = comparison) %>% 
    as_tibble(.)
  
}) %>% 
  set_names(comparison_list)

# read different comparisons - all genes
diffExp_all_results <- purrr::map(comparison_list, function(comparison){
  
  # read sheet
  read.xlsx(xlsxFile = diffExp_path_all, sheet = comparison) %>% 
    as_tibble(.)
  
}) %>% 
  set_names(comparison_list)

######################################################## MAIN CODE
### prepare data
## get protein coding genes
protein_coding_genes <- 
  genes_info_mouse %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id


### do the enrichment for all comparisons
go_enrich_list <- purrr::map(comparison_list, function(comparison){
  
  # get all results table for one comparison, filter protein coding genes
  diff_df_all <- 
    diffExp_all_results[[comparison]] %>% 
    left_join(., mouse_homologs, by = "gene_id") %>% 
    dplyr::filter(!is.na(mmusculus_homolog_ensembl_gene)) %>% 
    dplyr::filter(mmusculus_homolog_ensembl_gene %in% protein_coding_genes) %>% 
    dplyr::filter(!duplicated(mmusculus_homolog_ensembl_gene))
  
  # get significant results table for one comparison, filter protein coding genes
  diff_df_significant <- 
    diffExp_results[[comparison]] %>% 
    left_join(., mouse_homologs, by = "gene_id") %>% 
    dplyr::filter(!is.na(mmusculus_homolog_ensembl_gene)) %>% 
    dplyr::filter(mmusculus_homolog_ensembl_gene %in% protein_coding_genes)
  
  # prepare gene list
  geneList <- as.integer(diff_df_all$mmusculus_homolog_ensembl_gene %in% diff_df_significant$mmusculus_homolog_ensembl_gene)
  names(geneList) <- diff_df_all$mmusculus_homolog_ensembl_gene
  
  # make Probability Weighting Function for enrichment
  pwf_results <- nullp(DEgenes = geneList, genome = "mm10", id = "ensGene", plot.fit = FALSE)
  
  # do the enrichment
  goseq_results <- goseq(pwf = pwf_results, genome = "mm10", id = "ensGene")
  
  # filter results
  goseq_sign_results <- 
    goseq_results %>% 
    as_tibble(.) %>% 
    dplyr::filter(p.adjust(over_represented_pvalue, method = "BH") < 0.05)
  
  # save
  readr::write_csv(goseq_sign_results, 
                   file = file.path(outpath, str_c(experiment, "goseq.results", comparison, "mouse_homologs.significant.csv", sep = ".")))
  
  # return
  return(comparison)
  
})
