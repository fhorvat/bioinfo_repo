### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/phylogenetic_trees")

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

# animals table path
animals_tb_path <- file.path(inpath, "animal_list.20201012.csv")

# orthologs table path
orthologs_tb_path <- file.path(inpath, "mov10l1_orthologs.csv")

######################################################## READ DATA
# read annotation with protein IDs
animals_tb <- readr::read_csv(animals_tb_path)

# read orthologs table
orthologs_tb <- readr::read_csv(orthologs_tb_path)
  
######################################################## MAIN CODE
### prepare the data
# list of orthologs paper: 10.1093/gbe/evu105
# get scientific names
orthologs_tb %<>% 
  dplyr::mutate(scientific_name = str_remove_all(Species, ".*\\(|\\)"))

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 98, 97, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Aug 2020", "Sep 2019", "Jul 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://aug2020.archive.ensembl.org", 
                         "http://sep2019.archive.ensembl.org", 
                         "http://jul2019.archive.ensembl.org", 
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

# get animal names
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")

# list datasets
datasets_tb <- 
  listDatasets(mart) %>% 
  as_tibble(.)
  
# add dataset names to animals table
animals_tb_full <- 
  animals_tb %>% 
  dplyr::left_join(., datasets_tb, by = c("genome" = "version")) %>% 
  dplyr::left_join(., orthologs_tb, by = "scientific_name") %>% 
  dplyr::mutate(Type = str_remove(Type, "View Gene Tree"), 
                Orthologue = str_remove(Orthologue, "Compare Regions.*"), 
                gene_id = str_remove_all(Orthologue, ".*\\(|\\)"),
                gene_name = str_remove(Orthologue, "  \\(.*")) %>% 
  dplyr::select(common_name:dataset, 
                gene_id, gene_name,
                orthologue_type = Type, target_id = `Target %id`, query_id = `Query %id`, 
                coverage = `WGA Coverage`, confidence = `High Confidence`) %>% 
  dplyr::mutate(target_id = str_remove(target_id, " %") %>% as.numeric(.), 
                query_id = str_remove(target_id, " %") %>% as.numeric(.), 
                gene_id = replace(gene_id, scientific_name == "Homo sapiens", "ENSG00000073146"))

# get guys without annotated Mov10l1 gene
no_target <- 
  animals_tb_full %>% 
  dplyr::filter(is.na(gene_id))

# get guys with more than one orthologue
multiple_targets <-
  animals_tb_full %>% 
  dplyr::group_by(scientific_name) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::filter(query_id + target_id == max(target_id + query_id)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter((common_name != "dog" | gene_id == "ENSCAFG00030010927"))

# get guys with one orthologue 
one_target <- 
  animals_tb_full %>% 
  dplyr::filter(!(scientific_name %in% c(no_target$scientific_name, multiple_targets$scientific_name)))
  
# join all together
animals_tb_full <- rbind(one_target, multiple_targets, no_target)

# save
# readr::write_csv(animals_tb_full, file.path(outpath, "animal_list.ensembl_names.20201012.csv"))


### go through animals and download protein sequences from ENSEMBL
animals_gene_ids <- purrr::map(na.omit(animals_tb_full$dataset), function(ensembl_dataset){
  
  # load ENSEMBL mart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ensembl_dataset, host = "www.ensembl.org")
  # mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ensembl_dataset, host = "http://useast.ensembl.org")
  
  # set gene of interest
  gene_of_interest <- "ENSMUSG00000055839"
  
  # get gene homologs
  ensembl_homologs <-
    getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
          mart = mart) %>%
    as_tibble(.)
  
}) %>% 
  dplyr::bind_rows(.)


### go through animals and download protein sequences from ENSEMBL
animals_protein_seq <- purrr::map(na.omit(animals_tb_full$dataset), function(ensembl_dataset){

  ensembl_dataset <- "pmarinus_gene_ensembl"
  
  # load ENSEMBL mart
  # mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ensembl_dataset, host = ensembl_url)
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ensembl_dataset, host = "asia.ensembl.org")
  
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






