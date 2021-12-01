### INFO: 
### DATE: Tue May 28 15:30:02 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/KRAB_ZFP")

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
 
# set ensembl version and name
ensembl_version <- 93
ensembl_name <- "mmusculus"

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
genes_info %>% dplyr::filter(str_detect(gene_description, "ZFP"))

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

# get genes which have KRAB domain
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "pfam"), mart = mart) %>% 
  as_tibble(.) %>% 
  filter(pfam == "PF01352") %>% 
  dplyr::rename(gene_id = ensembl_gene_id) %>% 
  readr::write_csv(., path = file.path(outpath, str_c("ensembl", ensembl_version, "KRAB_domain.pfam01352.csv", sep = ".")))

