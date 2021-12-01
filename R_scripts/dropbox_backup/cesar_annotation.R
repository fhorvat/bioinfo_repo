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
# get mapping between human and mouse ENSEMBL gene id
mm_genes <- 
  getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
        filters = "ensembl_gene_id", values = genes_df$ENSMUSG, mart = mart_mm) %>% 
  dplyr::filter(hsapiens_homolog_ensembl_gene != "")

# get mapping between human ENSEMBL gene id and HGNC symbol
hs_genes <-
  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
        filters = "ensembl_gene_id", values = mm_genes$hsapiens_homolog_ensembl_gene, mart = mart_hs) %>% 
  dplyr::filter(hgnc_symbol != "")

# join tables
mm_hg_genes <- 
  inner_join(mm_genes, hs_genes, by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id")) %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarize(hgnc_symbol = str_c(unique(hgnc_symbol), collapse = ", "))

# upadate table
genes_df_updated <- 
  genes_df %>% 
  dplyr::left_join(., mm_hg_genes, by = c("ENSMUSG" = "ensembl_gene_id")) %>% 
  dplyr::select(1:3, hgnc_symbol_human_homolog = hgnc_symbol, everything())
 
# get missing values
genes_df_missing <- 
  genes_df_updated %>% 
  dplyr::filter(is.na(hgnc_symbol_human_homolog))

######################################################## 
# try another Ensembl version for missing values
# Ensembl versions
ens_ver <- 86
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

# get mapping between human and mouse ENSEMBL gene id
mm_genes <- 
  getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
        filters = "ensembl_gene_id", values = genes_df_missing$ENSMUSG, mart = mart_mm) %>% 
  dplyr::filter(hsapiens_homolog_ensembl_gene != "")

# get mapping between human ENSEMBL gene id and HGNC symbol
hs_genes <-
  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
        filters = "ensembl_gene_id", values = mm_genes$hsapiens_homolog_ensembl_gene, mart = mart_hs) %>% 
  dplyr::filter(hgnc_symbol != "")

# join tables
mm_hg_genes <- 
  inner_join(mm_genes, hs_genes, by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id")) %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarize(hgnc_symbol = str_c(unique(hgnc_symbol), collapse = ", "))

# upadate missing table
genes_df_updated_2 <- 
  genes_df_updated %>% 
  dplyr::left_join(., mm_hg_genes, by = c("ENSMUSG" = "ensembl_gene_id")) %>% 
  dplyr::mutate(hgnc_symbol_human_homolog = ifelse(is.na(hgnc_symbol_human_homolog), hgnc_symbol, hgnc_symbol_human_homolog))

# get missing values
genes_df_missing <- 
  genes_df_updated_2 %>% 
  dplyr::filter(is.na(hgnc_symbol_human_homolog))

# save
readr::write_csv(genes_df_updated_2, path = file.path(outpath, "180518_gene_list.ensembl.csv"))

# get mapping between mouse ENSEMBL gene id and HGNC symbol
getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
      filters = "ensembl_gene_id", values = "ENSMUSG00000061983", mart = mart_mm)

# get mapping between human ENSEMBL gene id and mouse homolog
getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "hgnc_symbol", values = "RPS12", mart = mart_hs)

# get mapping between human ENSEMBL gene id and mouse homolog
getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
      filters = "ensembl_gene_id", values = "ENSG00000112306", mart = mart_hs)

# MGI Symbol

# genes_df %>% 
#   dplyr::filter(ENSMUSG == "ENSMUSG00000024266")

# # get list of ensembl mouse strains
# mice <- 
#   useMart("ENSEMBL_MART_MOUSE", host = ensembl_url) %>% 
#   listDatasets(.) %>% 
#   dplyr::transmute(ensembl_id = dataset %>% str_remove(., "_gene_ensembl"),
#                    animal = description %>% str_remove(., " genes \\(.*"), 
#                    ensembl_genome = version)
# 
# # get list of ensembl animals
# animals <- 
#   useMart("ENSEMBL_MART_ENSEMBL", host = ensembl_url) %>% 
#   listDatasets(.) %>% 
#   dplyr::transmute(ensembl_id = dataset %>% str_remove(., "_gene_ensembl"),
#                    animal = description %>% str_remove(., " genes \\(.*"), 
#                    ensembl_genome = version) %>% 
#   rbind(., mice) %T>%
#   readr::write_csv(., path = file.path("/common/DB/genome_reference/scripts", "ensembl.91.genomes.csv"))
#
# # get mapping between human and mouse ENSEMBL gene id
# hg_mm_genes <- 
#   getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), mart = mart) %>% 
#   dplyr::filter(mmusculus_homolog_ensembl_gene %in% genes_df$ENSMUSG)
# 
# # get mapping between gene names and ensembl_gene_id
# hg_genes <- 
#   getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id",
#         values = hg_mm_genes$ensembl_gene_id, mart = mart) %>% 
#   dplyr::left_join(., hg_mm_genes, by = "ensembl_gene_id")

# # get info about genes
# ensembl_info <- 
#   getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"), 
#         filters = "ensembl_gene_id", values = ensembl_genes$ensembl_gene_id, mart = mart) %>% 
#   data.table::as.data.table(.) %>% 
#   data.table::setkey(., ensembl_gene_id)

# ### get info about homolog genes
# # filter animals
# animals_filt <- c("btaurus", "ocuniculus",
#                   "rnorvegicus", "ccrigri", "cchok1gshd", "sscrofa",
#                   "mmusculus")
# 
# # get info about homolog genes
# ids <- purrr::map(animals_filt, function(attr){
#   
#   homologs_df <- 
#     getBM(attributes = c("ensembl_gene_id", 
#                          str_c(attr, "_homolog_ensembl_gene"), 
#                          str_c(attr, "_homolog_associated_gene_name")),
#           filters = "ensembl_gene_id", values = ensembl_genes$ensembl_gene_id, mart = mart) %>% 
#     tibble::as.tibble(.) %>% 
#     dplyr::group_by(ensembl_gene_id) %>% 
#     dplyr::summarise_at(vars(contains("_homolog_associated_gene_name")), 
#                         funs(homolog_gene_name = str_c(unique(.), collapse = ", "))) %>% 
#     dplyr::filter(homolog_gene_name != "") %>% 
#     data.table::as.data.table(.) %>% 
#     data.table::setnames(., old = "homolog_gene_name", new = str_c(attr, ".homolog_gene_name")) %>% 
#     data.table::setkey(., ensembl_gene_id)
#   
#   return(homologs_df)
#   
# })
# 
# # join tables
# ensembl_hg <- 
#   purrr::reduce(.x = ids, .f = merge, by = "ensembl_gene_id") %>% 
#   merge(ensembl_info %>% dplyr::select(ensembl_gene_id, hgnc_symbol), .)
# 
# # filter by gene list
# ensembl_hg_filt <- 
#   ensembl_hg %>% 
#   dplyr::filter()
# 
# ensembl_info %>% 
#   dplyr::filter(hgnc_symbol %in% genes_list) %>% 
#   dplyr::arrange(hgnc_symbol)