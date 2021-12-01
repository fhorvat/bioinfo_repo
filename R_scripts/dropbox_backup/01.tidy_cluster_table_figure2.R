### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Shami_2020_DevCell_GSE142585")

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

library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 94

# table of consensus SPG clusters
clusters_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Shami_2020_DevCell_GSE142585/Analysis/expression"
clusters_tb_path <- list.files(clusters_path, "-mmc4\\.xlsx", full.names = T)

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_tb_path, sheet = 3) %>% 
  as_tibble(.)

######################################################## MAIN CODE
## load mart from Ensembl Biomart
### download gene info from Ensembl
# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(100, 99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Apr 2020", "Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://jan2020.archive.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive


### get gene IDs of human-mouse 1-1 orthologs
# set animal
ensembl_name <- "hsapiens"

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get gene names
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"), mart = mart) %>% 
  as_tibble(.) 

# get info about genes
ensembl_info_human <- 
  ensembl_info %>% 
  dplyr::filter(mmusculus_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::select(human_gene_name = external_gene_name, human_gene_id = ensembl_gene_id, mouse_gene_id = mmusculus_homolog_ensembl_gene)
  

### get gene IDs of mouse-hamster 1-1 orthologs
# set animal
ensembl_name <- "mmusculus"

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get gene names
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "mauratus_homolog_ensembl_gene", "mauratus_homolog_orthology_type"), mart = mart) %>% 
  as_tibble(.) 

# get info about genes, save
ensembl_info_mouse <- 
  ensembl_info %>% 
  dplyr::filter(mauratus_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::select(mouse_gene_id = ensembl_gene_id, hamster_gene_id = mauratus_homolog_ensembl_gene)


### join with original table
# cluster table
clusters_tb_hamster <- 
  clusters_tb %>% 
  dplyr::left_join(., ensembl_info_human, by = c("X1" = "human_gene_name")) %>% 
  dplyr::left_join(., ensembl_info_mouse, by = "mouse_gene_id") %>% 
  dplyr::rename(human_gene_name = X1) %T>%
  readr::write_csv(., file.path(outpath, "Shami_2020_DevCell_GSE142585.figure2_clusters.hamster_gene_id.csv"))
  