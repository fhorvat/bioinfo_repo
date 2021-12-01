### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Mollusca_project/orhtologs")

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

library(Biostrings)
library(taxize)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# annotation table path
annotation_path <- file.path(inpath, "Dicer1.ribonuclease_III_domain.Metazoa.OrthoDB.20211117.table.txt")

# protein sequences path
fasta_path <- file.path(inpath, "Dicer1.ribonuclease_III_domain.Metazoa.OrthoDB.20211117.fasta")

######################################################## READ DATA
# read annotation with protein IDs
annotation_tb <- readr::read_delim(annotation_path, delim = "\t")

# read protein sequences
protein_seq <- Biostrings::readAAStringSet(filepath = fasta_path)

######################################################## MAIN CODE
### for each animal get classification
# get tax_id
tax_id <- 
  str_remove(annotation_tb$organism_taxid, "_.*") %>% 
  unique(.)

# get classification list
class_list <- taxize::classification(tax_id, db = "ncbi")

# clean list, transform to table
class_tb <- purrr::map(names(class_list), function(tax_id){
  
  # get one table, set name, clean
  class_list[[tax_id]] %>% 
    as_tibble(.) %>% 
    dplyr::mutate(organism_taxid = tax_id) 
  
}) %>% 
  dplyr::bind_rows(.)

# transform to wide table
class_tb_wide <- 
  class_tb %>% 
  dplyr::filter(rank %in% c("phylum", "subphylum", "class")) %>% 
  dplyr::select(organism_taxid, name, rank) %>% 
  tidyr::pivot_wider(id_cols = organism_taxid, names_from = rank, values_from = name)

# join with annotation 
annotation_tb_class <- 
  annotation_tb %>% 
  dplyr::mutate(organism_taxid = str_remove(organism_taxid, "_.*")) %>% 
  dplyr::left_join(., class_tb_wide, by = "organism_taxid") %>% 
  dplyr::select(organism_name, phylum, subphylum, class, organism_taxid, int_prot_id, pub_gene_id) %>% 
  dplyr::mutate(phylum = factor(phylum, levels = c("Placozoa", "Porifera", "Cnidaria",
                                                   "Echinodermata", "Hemichordata", "Chordata", 
                                                   "Mollusca", "Annelida", "Brachiopoda", "Platyhelminthes", 
                                                   "Priapulida", "Nematoda", "Arthropoda"))) %>% 
  dplyr::arrange(phylum, subphylum, class, organism_name)

# clean sequence names
names(protein_seq) <- str_remove(names(protein_seq), " .*")

# create protein table
protein_tb <- tibble(int_prot_id = names(protein_seq), 
                     protein_seq = as.character(protein_seq))

# join with annotation table
protein_annotation_tb <- 
  annotation_tb_class %>% 
  dplyr::left_join(., protein_tb, by = "int_prot_id")

# save
openxlsx::write.xlsx(x = protein_annotation_tb, 
                     file = file.path(outpath, "Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences.xlsx"),
                     asTable = T)
