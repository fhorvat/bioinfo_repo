### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees/sequences")

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
library(GenomicRanges)
library(Biostrings)
library(msa)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# annotation with protein IDs path
annotation_proteins_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees/sequences/Piwi_Mov10l1.ensembl.proteins.csv"

######################################################## READ DATA
# read annotation with protein IDs
annotation_proteins <- readr::read_csv(annotation_proteins_path)

######################################################## MAIN CODE
# # annotation path
# annotation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_proteins_expression/results/piRNA_proteins.gene_ids.with_notes.csv"
# 
# # read annotation
# annotation_tb <- readr::read_csv(annotation_path)
# 
# # get PIWI proteins and Mov10l1
# annotation_tb_filt <- 
#   annotation_tb %>% 
#   dplyr::filter(str_detect(gene_name, "Piwil|Mov10l1")) %>% 
#   dplyr::select(-note) %>% 
#   tidyr::pivot_longer(., cols = -gene_name, names_prefix = "gene_id\\.", names_to = "animal", values_to = "gene_id") %>% 
#   dplyr::filter(!is.na(gene_id)) %>% 
#   dplyr::arrange(animal)
# 
# # save
# readr::write_csv(annotation_tb_filt, file.path(outpath, "Piwi_Mov10l1.ensembl.csv"))


### prepare data
# animals and their ensembl names
animals_tb <- tibble(animal = c("mouse", "rat", "golden_hamster", "cow", "human", "chicken"),
                     ensembl_name = c("mmusculus", "rnorvegicus", "mauratus", "btaurus", "hsapiens", "ggallus"))

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
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


### go through animals and download protein sequences from ENSEMBL
animals_protein_seq <- purrr::map(animals_tb$animal, function(animal_name){
  
  # filter ensembl name
  ensembl_name <- 
    animals_tb %>% 
    dplyr::filter(animal == animal_name) %$%
    ensembl_name
  
  # filter ensembl name
  annotation_animal <-
    annotation_proteins %>%
    dplyr::filter(animal == animal_name) %$%
    transcript_id
  
  # load ENSEMBL mart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)
  
  # get protein sequences
  protein_seq <- 
    getSequence(id = annotation_animal, type = "ensembl_transcript_id", 
                seqType = "peptide", mart = mart) %>% 
    as_tibble(.)
  
}) %>% 
  dplyr::bind_rows(.)

# add gene names
annotation_proteins_seq <- 
  annotation_proteins %>% 
  dplyr::left_join(., animals_protein_seq, by = c("transcript_id" = "ensembl_transcript_id")) %>% 
  # dplyr::filter(!is.na(peptide)) %>% 
  dplyr::mutate(protein_name = str_c(gene_name, animal, gene_id, transcript_id, sep = ".")) %>% 
  dplyr::select(protein_name, protein_seq = peptide, gene_name) %>% 
  dplyr::mutate(protein_seq = replace(protein_seq, is.na(protein_seq), "X"))

### filter and save
purrr::map(unique(annotation_proteins_seq$gene_name), function(gene){

  # filter
  protein_tb <- 
    annotation_proteins_seq %>% 
    dplyr::filter(gene_name == gene)
  
  # create AAStringSet
  protein_seq <- Biostrings::AAStringSet(protein_tb$protein_seq)
  names(protein_seq) <- protein_tb$protein_name 
  
  # save
  writeXStringSet(protein_seq, filepath = file.path(outpath, str_c(gene, ".AA.fasta")))

  # return
  return(gene)
  
})






