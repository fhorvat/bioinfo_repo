### INFO: 
### DATE: Mon Oct 22 13:09:08 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_expression/GO_terms")

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
library(clusterProfiler)
library(GOplot)
library(org.Mm.eg.db)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# path to differentialy expressed genes
diff_path <- file.path(inpath, "DE_genes.csv")

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
gene_info_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read list of diff. expressed genes
diff_df <- readr::read_csv(diff_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
### get GO terms from ENSEMBL
ensembl_release <- 93
ensembl_name <- "mmusculus"

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(93, 92, 91, 89, 86),
         date = c("Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_release) %$%
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
  getBM(attributes = c("ensembl_gene_id", "entrezgene", "go_id", "name_1006"), mart = mart) %>%
  as.tibble(.) %>%
  dplyr::rename(gene_id = ensembl_gene_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(entrezgene = str_c(unique(entrezgene), collapse = "/"),
                   go_id = str_c(unique(go_id), collapse = "/"),
                   go_name = str_c(unique(name_1006), collapse = "/")) %T>%
  readr::write_csv(., path = file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.GOterms.csv"))

# join with diff. expressed genes table
diff_df %<>%
  dplyr::mutate(gene_id = str_remove(gene_id, "\\..*$")) %>%
  dplyr::left_join(., ensembl_info, by = "gene_id") %T>%
  readr::write_csv(., file.path(outpath, "DE_genes.GO_terms.csv"))

