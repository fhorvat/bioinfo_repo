### INFO: CNOT gene annotation from Biomart
### DATE: Fri Mar 09 22:55:15 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(purrr)
library(ggplot2)

library(GenomicRanges)
library(biomaRt)

######################################################## PATH VARIABLES
outpath <- getwd()

# ensembl gene info path
geneInfo_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv" 

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read table with ensembl gene info 
geneInfo_df <- readr::read_csv(file = geneInfo_path)

######################################################## MAIN CODE
### get info about genes from Biomart
# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  count <- count + 1
  print(str_c("mmusculus_gene_ensembl", " ", count))
  mart <- tryCatch(expr = useMart(host = "may2017.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"),
                   error = function(x) return("error"))
}

# get info about genes
ensembl_gene_info <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", 
                                          "gene_biotype", "description"), mart = mart)

# get external info about genes
ensembl_ext_gene_info <- 
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "unigene", "ucsc"), mart = mart) %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarise(external_gene_name = str_c(unique(external_gene_name), collapse = ";"), 
                   refseq_mrna = str_c(unique(refseq_mrna), collapse = ";"), 
                   unigene = str_c(unique(unigene), collapse = ";"), 
                   ucsc = str_c(unique(ucsc), collapse = ";")) %>% 
  dplyr::mutate_at(.vars = 2:ncol(.), .funs = funs(str_replace_all(., ";{1,}$|^;{1,}", "")))

# get info about transcripts
ensembl_transcript_info <-
  getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), mart = mart) %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarise(ensembl_transcript_id = str_c(ensembl_transcript_id, collapse = ";")) %>% 
  dplyr::mutate_at(.vars = 2:ncol(.), .funs = funs(str_replace_all(., ";{1,}$|^;{1,}", "")))

# get info about probes
ensembl_probes_info <- 
  getBM(attributes = c("ensembl_gene_id", "affy_moe430a", "affy_moe430b"), mart = mart)  %>% 
  group_by(ensembl_gene_id) %>% 
  dplyr::summarise(affy_moe430a = str_c(affy_moe430a, collapse = ";"), 
                   affy_moe430b = str_c(affy_moe430b, collapse = ";")) %>% 
  dplyr::mutate_at(.vars = 2:ncol(.), .funs = funs(str_replace_all(., ";{1,}$|^;{1,}", "")))

# get info about homologs
mart <- "error"
count <- 0
while(class(mart) == "character"){
  count <- count + 1
  print(str_c("mmusculus_gene_ensembl", " ", count))
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"),
                   error = function(x) return("error"))
} 

ensembl_homologs_info <- getBM(attributes = c("ensembl_gene_id", "btaurus_homolog_ensembl_gene", "mauratus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"),
                               mart = mart)

### merge all, change, save
ensembl_gene_df <- 
  list(ensembl_gene_info, ensembl_ext_gene_info, ensembl_transcript_info, ensembl_homologs_info, ensembl_probes_info) %>% 
  purrr::reduce(., .f = left_join, by = "ensembl_gene_id") %>% 
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name, gene_description = description,
                seqnames = chromosome_name, start = start_position, end = end_position,
                btaurus_gene_id = btaurus_homolog_ensembl_gene, mauratus_gene_id = mauratus_homolog_ensembl_gene) %>%
  dplyr::mutate(seqnames = str_c("chr", seqnames),
                seqnames = replace(seqnames, seqnames == "chrMT", "chrM"),
                strand = replace(strand, strand == "1", "+"),
                strand = replace(strand, strand == "-1", "-"),
                strand = replace(strand, strand == "C", "*")) %T>%
  readr::write_csv(., path = file.path(outpath, "Mus_musculus.GRCm38.89.20180305.geneInfo.long.csv"))

# read
ensembl_gene_df <- readr::read_csv(file = file.path(outpath, "Mus_musculus.GRCm38.89.20180305.geneInfo.long.csv"))

# filter CNOT complex
CNOT_info <- 
  ensembl_gene_df %>% 
  dplyr::filter(str_detect(gene_name, "Cnot|Pan[2,3]{1}|Parn")) %>% 
  dplyr::mutate(gene_description = str_replace(gene_description, " \\[.*", "")) %>% 
  dplyr::mutate(coordinates = str_c(seqnames, ":", start, "-", end)) %>% 
  dplyr::select(gene_name, coordinates, strand, gene_id, everything()) %>% 
  dplyr::select(-c(seqnames:end, gene_description, gene_biotype)) %>% 
  dplyr::arrange(gene_name) %T>% 
  readr::write_csv(., path = file.path(outpath, "Mus_musculus.GRCm38.89.20180305.CNOTinfo.csv"))

  
