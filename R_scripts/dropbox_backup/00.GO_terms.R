### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
ensembl_version <- 99
gene_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
### get GO terms from ENSEMBL
ensembl_name <- "mmusculus"

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 93, 92, 91, 89, 86),
         date = c("Jan 2020", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://jan2020.archive.ensembl.org", 
                         "http://jul2018.archive.ensembl.org",
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
ensembl_info <-
  getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"), mart = mart) %>%
  as_tibble(.) %>%
  dplyr::rename(gene_id = ensembl_gene_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(go_id = str_c(unique(go_id), collapse = "/"),
                   go_name = str_c(unique(name_1006), collapse = "/")) %T>%
  readr::write_csv(., file = file.path(genome_dir, 
                                       gene_info_path %>% basename(.) %>% str_replace(., "geneInfo.csv", "GOterms.csv")))

