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
ensembl_version <- 99

# chosen genes path
genes_path <- list.files(inpath, "DevCell genes 200620.xlsx", full.names = T)

# mouse expression path
mouse_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Gan_2013_NatCommun_GSE35005/Analysis/expression"
mouse_tb_path <- list.files(mouse_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

# hamster expression path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression"
hamster_tb_path <- list.files(hamster_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

######################################################## READ DATA
# read genes table
genes_tb <- 
  openxlsx::read.xlsx(genes_path, sheet = 1) %>% 
  as_tibble(.)

# read mouse expression
mouse_tb <- readr::read_csv(mouse_tb_path)

# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

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


### get gene IDs of mouse-hamster 1-1 orthologs
# set animal
ensembl_name <- "mmusculus"

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get gene names
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mauratus_homolog_ensembl_gene", "mauratus_homolog_orthology_type"), mart = mart) %>% 
  as_tibble(.) 

# get info about genes, save
ensembl_info_mouse <- 
  ensembl_info %>% 
  # dplyr::filter(mauratus_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::select(mouse_gene_id = ensembl_gene_id, gene_name = external_gene_name, 
                hamster_gene_id = mauratus_homolog_ensembl_gene)


### join with original table
# genes table
genes_tb_ids <- 
  genes_tb %>% 
  dplyr::mutate(gene = str_to_sentence(gene)) %>% 
  dplyr::mutate(gene = replace(gene, gene == "Lin28", "Lin28a")) %>% 
  dplyr::select(gene_name = gene) %>% 
  dplyr::left_join(., ensembl_info_mouse, by = "gene_name") %>% 
  dplyr::mutate(hamster_gene_id = replace(hamster_gene_id, gene_name == "Piwil3", "gene14314")) %>% 
  dplyr::filter(hamster_gene_id != "ENSMAUG00000000336") %>% 
  dplyr::mutate(hamster_gene_id = replace(hamster_gene_id, hamster_gene_id == "", NA)) # remove other ortholog of FMR1

# join with mouse and hamster FPKM values
genes_tb_FPKM <- 
  genes_tb_ids %>% 
  dplyr::left_join(., mouse_tb %>% dplyr::select(mouse_gene_id = gene_id, primitive_SG_A, SG_A, SG_B, leptotene_SC, pachytene_SC, round_ST, elongative_ST), by = "mouse_gene_id") %>% 
  dplyr::left_join(., hamster_tb %>% dplyr::select(hamster_gene_id = gene_id, Mov10l_KO_13dpp:Mov10l_WT_adult), by = "hamster_gene_id") 

# save
readr::write_csv(genes_tb_FPKM, path = file.path(outpath, "DevCell_genes_200620.FPKM.Gan.hamster_testis.csv"))
