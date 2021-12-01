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

library(biomaRt)
library(goseq)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
experiment_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### results path
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", experiment)

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression.added_PIWIL3.stranded")
analysis_path <- file.path(analysis_path, "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq")

# differential expression results path
diffExp_path <- list.files(analysis_path, "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx", full.names = T)
diffExp_path_all <- list.files(analysis_path, "protein_coding.diffExp.DESeq2.genotype_age.all_results.xlsx", full.names = T)


### genome path
# set ensembl version
ensembl_version <- 99

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# genes info path
genes_info_path <- list.files(path = genome_path, 
                              pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), 
                              full.names = T)

# GO terms info
go_terms_info_path <- list.files(path = genome_path, 
                                 pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.GOterms.csv$"), 
                                 full.names = T)

######################################################## READ DATA
# read gene info
genes_info <- readr::read_csv(genes_info_path)

# read GO terms info
go_term_info <- readr::read_csv(go_terms_info_path)

### read list of diff. expressed genes
# set comparison
comparison_list <- c("Mov10l1_KO_8.5dpp_vs_Mov10l1_WT")

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
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id
  
# ## get transcripts info from biomaRt
# # load mart object
# mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mauratus_gene_ensembl", 
#                 host = "http://jan2020.archive.ensembl.org")
# 
# # get mean transcript length for each gene
# transcript_info <-
#   getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"), mart = mart) %>%
#   as_tibble(.) %>%
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::summarise(transcript_length = mean(transcript_length)) %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::rename(gene_id = ensembl_gene_id) %>% 
#   dplyr::filter(gene_id %in% protein_coding_genes)
# 
# # save
# readr::write_csv(transcript_info,
#                  file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.protein_coding.mean_transcript_length.csv"))

# read transcript length
transcript_info <- readr::read_csv(file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.protein_coding.mean_transcript_length.csv"))

## get gene to GO terms mappings
# filter table
go_term_info_tb <- 
  go_term_info %>% 
  dplyr::filter(gene_id %in% protein_coding_genes) 

# transform table to list
go_term_info_list <- str_split(go_term_info_tb$go_id, "/")
names(go_term_info_list) <- go_term_info_tb$gene_id


### do the enrichment for all comparisons
go_enrich_list <- purrr::map(comparison_list, function(comparison){
  
  # get all results table for one comparison, filter protein coding genes
  diff_df_all <- 
    diffExp_all_results[[comparison]] %>% 
    dplyr::filter(gene_id %in% protein_coding_genes) %>% 
    dplyr::left_join(., transcript_info, by = "gene_id")
  
  # get significant results table for one comparison, filter protein coding genes
  diff_df_significant <- 
    diffExp_results[[comparison]] %>% 
    dplyr::filter(gene_id %in% protein_coding_genes)

  # prepare gene list
  geneList <- as.integer(diff_df_all$gene_id %in% diff_df_significant$gene_id)
  names(geneList) <- diff_df_all$gene_id

  # make Probability Weighting Function for enrichment
  pwf_results <- nullp(geneList, bias.data = diff_df_all$transcript_length, plot.fit = FALSE)

  # do the enrichment
  goseq_results <- goseq(pwf = pwf_results, 
                         gene2cat = go_term_info_list)

  # filter results
  goseq_sign_results <- 
    goseq_results %>% 
    as_tibble(.) %>% 
    dplyr::filter(p.adjust(over_represented_pvalue, method = "BH") < 0.05)
  
  # save
  readr::write_csv(goseq_sign_results, 
                   file = file.path(outpath, str_c(experiment, "goseq.results", comparison, "significant.csv", sep = ".")))
  
  # return
  return(comparison)
  
})
