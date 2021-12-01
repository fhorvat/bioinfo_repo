### INFO: 
### DATE: Thu Nov 07 15:52:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/hackaton/prepare_data")

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


# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get hemoglobin gene ID (Hba-a1, hemoglobin alpha, adult chain 1)
hemoglobin_id <- "ENSG00000206172"

# homologs table
hemoglobin_homologs_tb <-
  listAttributes(mart) %>%
  as_tibble(.) %>%
  dplyr::filter(str_detect(name, "_homolog_ensembl_gene")) %>%
  dplyr::mutate(ensembl_name = str_remove(name, "_homolog_ensembl_gene")) %>%
  dplyr::filter(ensembl_name %in% animals_tb$ensembl_name)

# get gene IDs of homologs
homologs_tb <- purrr::map(hemoglobin_homologs_tb$ensembl_name, function(name){
  
  cat(name, "\n")
  
  # get info about genes
  ensembl_info <-
    getBM(attributes = c("ensembl_gene_id", str_c(name, c("_homolog_ensembl_gene", "_homolog_associated_gene_name"))),
          filters = "ensembl_gene_id",
          values = hemoglobin_id,
          mart = mart) %>%
    as_tibble(.) %>%
    magrittr::set_colnames(., str_remove(colnames(.), str_c(name, "_"))) %>%
    dplyr::mutate(ensembl_name = name)
  
  # return
  return(ensembl_info)
  
}) %>% 
  bind_rows(.) %>%
  dplyr::filter(!is.na(homolog_ensembl_gene))


### get sequences of genes and exon coordinates for each animal
ensembl_list <- purrr::map(unique(homologs_tb$ensembl_name), function(name){
  
  cat(name, "\n")
  
  # load ENSEMBL mart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(name, "_gene_ensembl"), host = ensembl_url)
  
  # filter table
  animal_info_tb <- homologs_tb %>% dplyr::filter(ensembl_name == name)
  
  # get gene sequence
  ensembl_seq <-
    getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                         "chromosome_name", "start_position", "end_position", "strand", "gene_exon_intron"),
          filters = "ensembl_gene_id",
          values = animal_info_tb$homolog_ensembl_gene,
          mart = mart) %>%
    as_tibble(.) %>%
    dplyr::mutate(ensembl_name = name, 
                  strand = ifelse(strand == 1, "+", "-"))
  
  ### get exon coordinates
  # getBM sometimes doesn't work and it returns an empty table
  ensembl_exons <- tibble()
  count <- 0
  
  while(nrow(ensembl_exons) == 0){
    
    count <- count + 1
    cat(str_c(name, " ", count), "\n")
    
    ensembl_exons <- 
      getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id", 
                           "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", 
                           "rank"),
            filters = "ensembl_gene_id",
            values = animal_info_tb$homolog_ensembl_gene,
            mart = mart) %>%
      as_tibble(.) %>%
      dplyr::arrange(ensembl_gene_id, rank) %>% 
      dplyr::mutate(ensembl_name = name, 
                    strand = ifelse(strand == 1, "+", "-"))
    
    # stop after 50 tries
    if(count == 50){
      stop("Something wrong with biomaRt!")
    }
    
  }

  
  # return list
  return(list(gene_seq = ensembl_seq, exon_coords = ensembl_exons))
  
}) %>% 
  set_names(., unique(homologs_tb$ensembl_name))


#### extract sequences and coordinates of the genes
### extract gene info
genes_list <- lapply(ensembl_list, "[[", 1)

### get all gene sequences on + strand (sequences are sense, so if the gene is on the - strand they need to be revComp)
# extract sequences
genes_seq <- 
  lapply(names(genes_list), function(name){
    
    # extract table
    genes_tb <- genes_list[[name]]
    
    # create an DNAStringSet
    genes_dna <- DNAStringSet(genes_tb$gene_exon_intron)
    
    # set names
    names(genes_dna) <- str_c(name, genes_tb$ensembl_gene_id, sep = ".")
    
    # return
    return(genes_dna)
    
  }) %>% 
  do.call(c, .)

# check
if(length(genes_seq) != nrow(homologs_tb)){
  stop("Something wrong with genes!")
}
 
# get indices of genes on minus strand 
minus_strands <- 
  lapply(genes_list, "[[", "strand") %>% 
  unlist(.) %>% 
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

# check
if(nrow(gene_tb) != nrow(homologs_tb)){
  stop("Something wrong with gene coordinates!")
}

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

# check
if(length(unique(exons_tidy$gene_id)) != nrow(homologs_tb)){
  stop("Something wrong with exonic coordinates!")
}
  

# get genes which have annotated alternatively spliced genes (exons repeating more than once)
alt_spliced_genes <- 
  exons_tidy %>%
  dplyr::group_by(ensembl_name, gene_id) %>%
  dplyr::filter(length(exon_rank) != length(unique(exon_rank))) %T>%
  write_csv(., file.path(outpath, "alt_spliced_genes.all_homologs.csv"))

# read table
constitutive_exons <- 
  lapply(list.files(inpath, "ExonsSpreadsheet*"), function(path){
    
    # read table
    exons <- 
      readr::read_csv(path) %>% 
      dplyr::filter(!is.na(No.)) %>% 
      dplyr::select(exon_id = `Exon / Intron`)
    
  }) %>% 
  dplyr::bind_rows(.)

# filter table
alt_spliced_genes_constitutive <- 
  alt_spliced_genes %>% 
  dplyr::filter(exon_id %in% constitutive_exons$exon_id)

# filter table
exons_final <-
  exons_tidy %>% 
  dplyr::filter(!(ensembl_name %in% alt_spliced_genes_constitutive$ensembl_name) | (exon_id %in% alt_spliced_genes_constitutive$exon_id)) %>% 
  dplyr::select(-c(exon_rank, exon_id))

###
# write data
readr::write_csv(exons_final, file.path(outpath, "hemoglobin.exonic_coordinates.all.csv"))
readr::write_csv(gene_tb, file.path(outpath, "hemoglobin.genomic_coordinates.all.csv"))
saveRDS(genes_seq, file.path(outpath, "hemoglobin.genomic_sequences.DNAStringSet.all.RDS"))

