### INFO: 
### DATE: Mon May 31 10:08:37 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Fugaku_FPKM")

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

# fpkm tables path
fpkm_path <- list.files(inpath, ".*\\.all\\.FPKM\\.csv", full.names = T)

######################################################## READ DATA
# read fpkm tables  
fpkm_list <- purrr::map(fpkm_path, function(path){
  
  # read
  readr::read_csv(file = path)
  
}) %>% 
  set_names(basename(fpkm_path))

######################################################## MAIN CODE
# set ensembl version and name
ensembl_version <- 93
ensembl_name <- "mmusculus"

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

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = "uswest.ensembl.org")

# get all attributes
ensembl_attr <- 
  biomaRt::listAttributes(mart = mart) %>% 
  as_tibble(.) %>% 
  dplyr::filter(page == "feature_page")

# get info about genes
ensembl_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)


