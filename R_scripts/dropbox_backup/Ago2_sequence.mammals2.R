### INFO: get sequence of Ago2 from different rodents
### DATE: 06. 10. 2017.
### AUTHOR: Vedran Franke, Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/Sequences")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(biomaRt)
library(Biostrings)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# set animals from which you want sequences 
animals_list <- c("guinea pig" = "cporcellus",
                  "rat" = "rnorvegicus",
                  "mouse" = "mmusculus",
                  "rabbit" = "ocuniculus",
                  # "squirrle" = "itridecemlineatus",
                  "kangaroo rat" = "dordii",
                  "human" = "hsapiens",
                  "pig" = "sscrofa"
                  # "elephant" = "lafricana",
                  # "bat" = "mlucifugus",
                  # "platypus" = "oanatinus",
                  # "dog" = "cfamiliaris",
                  # "cow" = "btaurus"
)

# list of animal assemblies
animals <- 
  useMart("ensembl") %>% 
  listDatasets(.) %$% 
  dataset %>% 
  str_replace(., "_gene_ensembl", "")

# filter animals by list
animals <- animals_list[animals_list %in% animals]

######################################################## MAIN CODE
### get Ensembl gene IDs of homolog of chosen gene in all mouse strains from Biomart
# set geneID
gene_symbol <- "Ago2"

# load Mart of any mouse strain database from ENSEMBL_MART_MOUSE (here I use 129S1/SvImJ = m129s1svimj)
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  count <- count + 1
  print(str_c("hsapiens_gene_ensembl", " ", count))
  mart <- tryCatch(expr = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                   error = function(x) return("error"))
}

# set attributes to fetch from loaded Mart
attrs <- 
  str_c(animals, "_homolog_ensembl_gene") %>% 
  c(., "ensembl_gene_id") %>% 
  .[. != "hsapiens_homolog_ensembl_gene"]

# download attributes
# getBM runs in loop because it is possible to fetch only 6 attributes in one call
set <- rep(1:((length(attrs) %/% 3) + 1), each = 3)[1:length(attrs)]
ids <- 
  lapply(unique(set), function(sub_set){
    homolog_geneID <- getBM(attributes = c(attrs[sub_set == set]), filters = c("external_gene_name"), values = gene_symbol, mart = mart)  
    return(homolog_geneID)
  }) %>% 
  unlist(.) %>% 
  magrittr::set_names(., ifelse(names(.) == "ensembl_gene_id", "hsapiens", str_replace(names(.), "_.+", ""))) %>% 
  .[!duplicated(.)]


### get sequences and coordinates of homolog of chosen gene from other species
# set datasets to load from Biomart
ensembl_datasets <- str_c(names(ids), "_gene_ensembl", sep = "")

# loops over the organisms and gets coordinates/sequences
bm_data <- 
  lapply(ensembl_datasets, function(ensembl_dataset){
    
    # load Mart
    mart <- "error"
    count <- 0
    while(class(mart) == "character"){
      count <- count + 1
      print(str_c(ensembl_dataset, " ", count))
      mart <- tryCatch(expr = useMart(biomart =  "ensembl", dataset  = ensembl_dataset),
                       error = function(x) return("error"))
    }
    
    # get cooridnates of exons
    exons_coords <- 
      getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_transcript_id", "strand"),
            filters = c("ensembl_gene_id"),
            values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")],
            mart = mart) %>% 
      dplyr::arrange(exon_chrom_start) %>% 
      dplyr::mutate(dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = ""))
    
    # get coordinates and whole sequences of genes
    genes_coords_seq <- 
      getBM(attributes = c("start_position", "end_position", "gene_exon_intron"),
            filters = c("ensembl_gene_id"),
            values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")],
            mart = mart) %>% 
      dplyr::mutate(dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = ""))
    
    return(list(exons_coords, genes_coords_seq))
    
  })

biomaRt::getSequence

### get sequences of chosen intron in gene from all mouse strains
# set intron in form LeftExonNumber_RightExonNumber (in sense direction)
intron_filt <- "1_2"

# get coordinates of introns from coordinates of exons
introns_coord <- 
  lapply(bm_data, function(x) x[[1]]) %>% 
  dplyr::bind_rows(.) %>% 
  group_by(dataset, ensembl_transcript_id) %>% 
  dplyr::select(chromosome_name, exon_chrom_start, exon_chrom_end, strand, dataset, ensembl_transcript_id) %>%
  dplyr::arrange(ensembl_transcript_id, exon_chrom_start) %>% 
  dplyr::mutate(start_int = exon_chrom_end + 1,
                end_int = lead(exon_chrom_start) - 1, 
                exon_index = ifelse(strand == "-1", n():1, 1:n())) %>%  
  dplyr::filter(!is.na(end_int)) %>%
  dplyr::mutate(intron_index = ifelse(strand == -1, exon_index - 1, exon_index + 1),
                intron_index = ifelse(strand == -1, str_c(intron_index, "_", exon_index), str_c(exon_index, "_", intron_index))) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-exon_chrom_start, -exon_chrom_end, -exon_index, -ensembl_transcript_id) %>% 
  dplyr::select(chromosome_name, start = start_int, end = end_int, everything()) %>% 
  dplyr::distinct(.) %>% 
  dplyr::left_join(., tibble(dataset = animals,
                             animal_name = names(animals)), 
                   by = "dataset") %>% 
  dplyr::filter(intron_index == intron_filt)

# get coordinates and sequences of genes
genes_coords_seq <- 
  lapply(bm_data, function(x) x[[2]]) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(gene_start = start_position, gene_end = end_position, gene_exon_intron, dataset)

# join with intron coordinates with gene coordinates, get coordinates, subset gene sequences to chosen intron sequences
seq_data <- 
  dplyr::left_join(introns_coord, genes_coords_seq, by = "dataset") %>% 
  dplyr::mutate(subset_start = ifelse(strand == "-1", gene_end - end + 1, start - gene_start + 2), 
                subset_end = subset_start + (end - start - 1), 
                gene_exon_intron = stringr::str_sub(string = gene_exon_intron, start = subset_start, end = subset_end)) %>% 
  dplyr::select(-c(gene_start, gene_end, subset_start, subset_end)) 


### write sequences
# create dirs
outdir <- file.path(outpath, gene_symbol)
dir.create(outdir, showWarnings = F)

# get sequences
ensembl_seq <- 
  seq_data %>% 
  dplyr::filter(gene_exon_intron != "Sequence unavailable") %$% 
  gene_exon_intron %>% 
  DNAStringSet(.) 

# set names, order
names(ensembl_seq) <- 
  str_replace(seq_data$dataset, "(_|\\.).+", "") %>% 
  str_c(., gene_symbol, str_c("intron", seq_data$intron_index), sep = '_')

# write sequences
ensembl_seq %<>% 
  .[order(width(.))] %>% 
  .[(letterFrequency(., letters = "N", as.prob = T) %>% as.vector(.)) < 1] %T>%
  writeXStringSet(., file.path(outdir, str_c(gene_symbol, "mouse_strains", str_c("intron", intron_filt), length(.), "fa", sep = ".")), format = "fasta")

