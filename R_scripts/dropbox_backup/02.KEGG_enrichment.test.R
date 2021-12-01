### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression/GO_terms_enrichment")

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
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(biomaRt)

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
ensembl_name <- "mmusculus"
ensembl_version <- 91
gene_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read gene info
gene_info <- 
  readr::read_csv(gene_info_path) %>% 
  dplyr::filter(gene_biotype == "protein_coding")

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

######################################################## MAIN CODE
### get additional info about genes from Biomart
## Ensembl versions
ensembl_url <-
  tibble(ens_version = c(93, 92, 91, 89, 86),
         date = c("Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://jul2018.archive.ensembl.org", 
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive

# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  
  count <- count + 1
  print(str_c(ensembl_name, "_gene_ensembl", " ", count))
  
  # load ENSEMBL mart
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
                   error = function(x) return("error"))
  
  # if error try mouse strains mart
  if(class(mart) == "character"){
    mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_MOUSE", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
                     error = function(x) return("error"))
  }
  
  # stop if count get too big
  if(count > 2){
    stop("Something's not right")
  }
  
}

# get info about genes
entrezid_info <- 
  getBM(attributes = c("ensembl_gene_id", "uniprot_gn"), mart = mart) %>% 
  as.tibble(.) %>% 
  dplyr::rename(gene_id = ensembl_gene_id)


### do the enrichment for both stages
go_enrich_list <- purrr::map(stages, function(stage){
  
  stage <- "GV"
  
  # get results table for one stage
  diff_df <- 
    diffExp_results[[stage]] %>% 
    dplyr::left_join(., entrezid_info, by = "gene_id")
  
  # KEGG enrichment
  kk <- enrichKEGG(gene = entrezid_info$uniprot_gn,
                   keyType = "uniprot", 
                   organism = "msu",
                   pvalueCutoff = 0.05)
  
  
}) %>%
  set_names(., stages)

