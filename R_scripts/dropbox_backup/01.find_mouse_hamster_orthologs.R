### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse")

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

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
## load mart from Ensembl Biomart
# set animal
ensembl_name <- "mmusculus"

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

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get info about genes, save
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "mauratus_homolog_ensembl_gene", "mauratus_homolog_orthology_type"), mart = mart) %>% 
  as_tibble(.) %>% 
  dplyr::filter(mauratus_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(ensembl_gene_id = gene_id, gene_name), by = "ensembl_gene_id") %>% 
  dplyr::select(hamster_gene_id = mauratus_homolog_ensembl_gene, 
                mouse_gene_id = ensembl_gene_id, 
                mouse_gene_name = gene_name) %T>% 
  readr::write_csv(., file.path(outpath, str_c("ensembl", ensembl_version, "mouse_vs_goldHamster.one2one_homologs.csv", sep = ".")))

