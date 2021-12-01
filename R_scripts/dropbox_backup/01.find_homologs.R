### INFO: 
### DATE: Thu Nov 07 15:52:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Elob_MSA/sequences")

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
library(Biostrings)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

######################################################## READ DATA

######################################################## MAIN CODE
### get additional info about genes from Biomart
## Ensembl versions
ensembl_url <-
  tibble(ens_version = c(98, 97, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Sep 2019", "Jul 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://sep2019.archive.ensembl.org", 
                         "http://jul2019.archive.ensembl.org", 
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

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = ensembl_url)

# homologs table
homologs_tb <- 
  listAttributes(mart) %>% 
  as_tibble(.) %>% 
  dplyr::filter(str_detect(name, "_homolog_ensembl_gene")) %>% 
  dplyr::select(name) %>% 
  dplyr::mutate(name = str_remove(name, "_homolog_ensembl_gene"))

## find homologs in other species
# set gene of interest
gene_of_interest <- "ENSMUSG00000055839"

# get gene homologs
ensembl_homologs <-
  getBM(attributes = c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene",
                       "mauratus_homolog_ensembl_gene", "btaurus_homolog_ensembl_gene"),
        filters = "ensembl_gene_id",
        values = gene_of_interest,
        mart = mart) %>%
  as_tibble(.)

# get transcripts




