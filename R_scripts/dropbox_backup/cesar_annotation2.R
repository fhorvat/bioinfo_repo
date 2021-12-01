### INFO: 
### DATE: Thu May 17 16:33:03 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/CESAR_annotation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)

library(biomaRt)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
simpleCap <- function(x) {
  paste(toupper(substring(x, 1, 1)), tolower(substring(x, 2)), sep = "")
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# annotation path
annot_path <- file.path(inpath, "approved_gene_names.xlsx")

# genes path
genes_path <- file.path(inpath, "180517_gene_list.xlsx")

######################################################## READ DATA
# annotation
annot_df <- readxl::read_xlsx(path = annot_path, col_names = "hg.gene_name")

# list of genes
genes_df <- 
  readxl::read_xlsx(path = genes_path) %>%
  dplyr::mutate_all(funs(stringr::str_trim(., side = "both")))

# Ensembl versions
ens_ver <- 92
ensembl_url <- 
  tibble(ens_version = c(92, 91, 86), 
         date = c("Apr 2018", "91 Dec 2017", "Oct 2016"), 
         URL_archive = c("http://apr2018.archive.ensembl.org", 
                         "http://dec2017.archive.ensembl.org", 
                         "http://oct2016.archive.ensembl.org")) %>% 
  dplyr::filter(ens_version == ens_ver) %$% 
  URL_archive

# load Mart of human database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
ensembl_species <- "mmusculus"
mart_mm <- "error"
count <- 0
while(class(mart_mm) == "character"){
  count <- count + 1
  print(str_c(ensembl_species, "_gene_ensembl", " ", count))
  mart_mm <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = str_c(ensembl_species, "_gene_ensembl"), 
                                     host = ensembl_url), 
                      error = function(x) return("error"))
}

# load Mart of human database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
ensembl_species <- "hsapiens"
mart_hs <- "error"
count <- 0
while(class(mart_hs) == "character"){
  count <- count + 1
  print(str_c(ensembl_species, "_gene_ensembl", " ", count))
  mart_hs <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = str_c(ensembl_species, "_gene_ensembl"),
                                     host = ensembl_url),
                      error = function(x) return("error"))
}

######################################################## MAIN CODE
# clean mgi_symbol
mgi_df <- 
  genes_df %>% 
  dplyr::select(ID, official_name = `official name`) %>% 
  tidyr::separate(col = official_name, into = c("official_name1", "official_name2"), sep = "/") %>% 
  tidyr::gather(name, official_name, -ID) %>% 
  dplyr::select(-name) %>% 
  dplyr::filter(!is.na(official_name)) %>% 
  dplyr::mutate(official_name = simpleCap(official_name)) %>% 
  dplyr::select(ID, official_name_clean = official_name)

# get mapping between MGI symbol and ensembly gene id
mm_genes <- 
  getBM(attributes = c("mgi_symbol", "ensembl_gene_id"), 
        filters = "mgi_symbol", values = mgi_df$official_name_clean, mart = mart_mm) %>% 
  dplyr::left_join(., mgi_df, by = c("mgi_symbol" = "official_name_clean"))

# get human homologs
mm_hs_genes <- 
  getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
        filters = "ensembl_gene_id", values = mm_genes$ensembl_gene_id, mart = mart_mm) %>% 
  dplyr::filter(hsapiens_homolog_ensembl_gene != "")

# get HGNC_symbol of human homolog
hs_genes <- 
  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
        filters = "ensembl_gene_id", values = mm_hs_genes$hsapiens_homolog_ensembl_gene, mart = mart_hs) %>% 
  dplyr::filter(hgnc_symbol != "")

### join
mm_hs_genes %<>%
  dplyr::left_join(., hs_genes, by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id"))

mm_genes %<>%
  dplyr::left_join(., mm_hs_genes, by = "ensembl_gene_id")

mgi_df %<>%
  dplyr::left_join(., mm_genes, by = "ID") 

### clean and save
# filter by whol Hiller data set (keep only those which are present in Hiller)
mgi_df_clean <- 
  mgi_df %>% 
  dplyr::filter(hgnc_symbol %in% annot_df$hg.gene_name)

# add to original table
genes_df_clean <- 
  genes_df %>% 
  dplyr::left_join(., mgi_df_clean, by = "ID") %>% 
  dplyr::filter(!duplicated(ID)) %>% 
  dplyr::select(unique_ID = ID, hgnc_symbol, `function`, `official name`, mgi_symbol, `synonyms (from biogps)`, ENSMUSG, refseq, note, 
                ensembl_gene_id, hsapiens_homolog_ensembl_gene)

# save
readr::write_csv(x = genes_df_clean, path = file.path(inpath, "180517_gene_list.HGNC_symbol.csv"))
  