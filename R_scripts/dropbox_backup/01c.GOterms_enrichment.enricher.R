### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(openxlsx)
library(enrichR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# path to results of differential expression analysis
diffExp_results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression/results/Lnc1_KO.2018_Dec.diffExp.GRCm38.91.protein_coding.significant.xlsx"

# genome_dir path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
ensembl_version <- 91
gene_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# GO terms info
go_terms_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.GOterms.csv$"), full.names = T)

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read GO terms info
go_term_info <- readr::read_csv(go_terms_info_path)

### read list of diff. expressed genes
# set stages
stages <- c("GV", "MII")

# read different stages
diffExp_results <- purrr::map(stages, function(stage){
  
  # read sheet
  read.xlsx(xlsxFile = diffExp_results_path, sheet = stage) %>% 
    as.tibble(.)
  
}) %>% 
  set_names(stages)

# list data in enricher database
# dbs <- listEnrichrDbs()

######################################################## MAIN CODE
### do the enrichment for both stages
go_enrich_list <- purrr::map(stages, function(stage){
  
  stage <- "GV"
  
  # get results table for one stage
  diff_df <- diffExp_results[[stage]]
  
  dbs_chosen <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Mouse")
  enriched_list <- enrichr(diff_df$gene_name, dbs_chosen) 
  
})
