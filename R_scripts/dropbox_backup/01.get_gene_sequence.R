### INFO: 
### DATE: Thu Nov 07 15:52:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test")

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
library(Biostrings)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 97
ensembl_name <- "hsapiens"

# list of animals path
animals_path <- file.path(inpath, "animal_list.csv")

######################################################## READ DATA
# read animals table
animals_tb <- readr::read_csv(animals_path)

######################################################## MAIN CODE
### get additional info about genes from Biomart
## Ensembl versions
ensembl_url <-
  tibble(ens_version = c(98, 97, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Sep 2019", "Jul 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://sep2019.archive.ensembl.org", 
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


# # load ENSEMBL mart
# mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)
# 
# # get hemoglobin gene ID (Hba-a1, hemoglobin alpha, adult chain 1)
# hemoglobin_id <- "ENSG00000206172"
# 
# # homologs table
# homologs_tb <-
#   listAttributes(mart) %>%
#   as_tibble(.) %>%
#   dplyr::filter(str_detect(name, "_homolog_ensembl_gene")) %>%
#   dplyr::mutate(animal = str_remove(name, "_homolog_ensembl_gene")) %>%
#   dplyr::filter(animal %in% animals_tb$ensembl_name) %>% 
#   dplyr::filter(animal %in% c("hmale", "hfemale"))
# 
# # get gene IDs of homologs
# homologs_list <- purrr::map(homologs_tb$animal, function(animal){
# 
#   # get info about genes
#   ensembl_info <-
#     getBM(attributes = c("ensembl_gene_id", str_c(animal, c("_homolog_ensembl_gene", "_homolog_associated_gene_name"))),
#           filters = "ensembl_gene_id",
#           values = hemoglobin_id,
#           mart = mart) %>%
#     as_tibble(.) %>%
#     magrittr::set_colnames(., str_remove(colnames(.), str_c(animal, "_"))) %>%
#     dplyr::mutate(animal = animal)
# 
#   # return
#   return(ensembl_info)
# 
# })
# 
# # create table
# homologs_tb <-
#   homologs_list %>%
#   bind_rows(.) %>%
#   dplyr::filter(!is.na(homolog_ensembl_gene)) %>%
#   dplyr::left_join(., animals_tb, by = c("animal" = "ensembl_name"))


### get sequences of genes and exon coordinates for each animal
ensembl_list <- purrr::map(animals_tb$ensembl_name, function(ensembl_name){
  
  cat(ensembl_name, "\n")
  
  # load ENSEMBL mart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)
  
  # filter table
  animal_info_tb <- animals_tb %>% dplyr::filter(ensembl_name == ensembl_name)
  
  # get gene sequence
  ensembl_seq <-
    getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                         "chromosome_name", "start_position", "end_position", "strand", "gene_exon_intron"),
          filters = "ensembl_gene_id",
          values = animal_info_tb$ensembl_gene_id,
          mart = mart) %>%
    as_tibble(.) %>%
    dplyr::mutate(ensembl_name = ensembl_name, 
                  strand = ifelse(strand == 1, "+", "-"))
  
  # get exon coordinates
  ensembl_exons <- 
    getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id", 
                         "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", 
                         "rank",  "gene_exon"),
          filters = "ensembl_gene_id",
          values = animal_info_tb$ensembl_gene_id,
          mart = mart) %>%
    as_tibble(.) %>%
    dplyr::mutate(ensembl_name = ensembl_name, 
                  strand = ifelse(strand == 1, "+", "-"))
  
  # return list
  return(list(gene_seq = ensembl_seq, exon_coords = ensembl_exons))
  
})

# set names
names(ensembl_list) <- animals_tb$ensembl_name


#### extract sequences and coordinates of the genes
### extract gene info
genes_list <- lapply(ensembl_list, "[[", 1)

### get all gene sequences on + strand (sequences are sense, so if the gene is on the - strand they need to be revComp)
# extract sequences
genes_seq <- 
  lapply(genes_list, function(x) DNAString(x$gene_exon_intron)) %>% 
  DNAStringSet(.)

# get indices of genes on minus strand 
minus_strands <- 
  sapply(genes_list, "[[", "strand") %>% 
  magrittr::equals(., "-") %>% 
  which(.)

# reverse complement
genes_seq[minus_strands] <- reverseComplement(genes_seq[minus_strands])

## get all gene coordinates
gene_tb <- 
  lapply(genes_list, function(gene_tb){
    
    # tidy table
    gene_tb %>% 
      dplyr::select(seqnames = chromosome_name,
                    start = start_position, 
                    end = end_position, 
                    gene_id = ensembl_gene_id, 
                    ensembl_name) %>% 
      dplyr::mutate_if(is.numeric, as.character)
    
  }) %>% 
  bind_rows(.)


### get all constitutive exons
# extract exons
exons_list <- lapply(ensembl_list, "[[", 2)

# extract tidy exon coordinates
exons_tidy <- 
  lapply(exons_list, function(exon){
    
    # tidy table
    exon %>% 
      dplyr::select(seqnames = chromosome_name, start = exon_chrom_start, end = exon_chrom_end, 
                    gene_id = ensembl_gene_id, exon_id = ensembl_exon_id, exon_rank = rank, 
                    ensembl_name) %>% 
      dplyr::mutate_if(is.numeric, as.character)
    
  }) %>% 
  bind_rows(.) 

# # write table with genes which have annotated alternatively spliced genes (exons repeating more than once)
# exons_tidy %>% 
#   dplyr::group_by(ensembl_name) %>% 
#   dplyr::filter(length(exon_rank) != length(unique(exon_rank))) %>% 
#   write_csv(., file.path(outpath, "multiple_exons.csv"))

# read table
constitutive_exons <- readr::read_csv(file.path(inpath, "constitutive_exons.csv"))

# filter table
exons_tidy %<>% 
  dplyr::filter(!(ensembl_name %in% constitutive_exons$ensembl_name) | (exon_id %in% constitutive_exons$exon_id)) %>% 
  dplyr::select(-c(exon_rank, exon_id))


###
# write data
readr::write_csv(exons_tidy, file.path(outpath, "hemoglobin_exons_coordinates.csv"))
readr::write_csv(gene_tb, file.path(outpath, "hemoglobin_gene_coordinates.csv"))
saveRDS(genes_seq, file.path(outpath, "hemoglobin_gene_sequences.DNAStringSet.RDS"))

