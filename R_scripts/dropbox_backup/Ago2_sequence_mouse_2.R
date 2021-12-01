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

library(data.table)
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
# list of mouse strains
animals <- 
  useMart("ENSEMBL_MART_MOUSE") %>% 
  listDatasets(.)

animals <- 
  animals %$% 
  dataset %>% 
  str_replace(., "_gene_ensembl", "") %>% 
  set_names(., animals %$%
              description %>% 
              str_replace_all(., "Mouse | genes.*", "")) %>% 
  c(., "mmusculus" %>% set_names("Mus_musculus"))

######################################################## MAIN CODE
# set geneID
gene_symbol <- "Ago2"
gene_id <- "MGP_129S1SvImJ_G0021919"

### get geneIDs of homolog ensembl genes from other species
attrs <- 
  str_c(animals, "_homolog_ensembl_gene") %>% 
  .[. != "m129s1svimj_homolog_ensembl_gene"]

## sometimes function useMart isn't able to connect to server and returns error
# this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  count <- count + 1
  print(str_c("m129s1svimj_gene_ensembl", " ", count))
  mart <- tryCatch(expr = useMart("ENSEMBL_MART_MOUSE", dataset = "m129s1svimj_gene_ensembl"), 
                   error = function(x) return("error"))
}

## because it is possible to fetch only 6 attributes in one call, getBM runs in loop
set <- rep(1:((length(attrs) %/% 3) + 1), each = 3)[1:length(attrs)]
ids <- 
  lapply(unique(set), function(sub_set){
    homolog_geneID <- getBM(attributes = c(attrs[sub_set == set]), filters = c("ensembl_gene_id"), values = gene_id, mart = mart)  
    return(homolog_geneID)
}) %>% 
  unlist(.) %>% 
  magrittr::set_names(str_replace(names(.),'_.+','')) %>% 
  .[!duplicated(.)]


### get sequences of homolog ensembl genes from other species
ensembl_datasets <- str_c(animals, "_gene_ensembl", sep = "")

## loops over the organisms and fetches sequences
ensembl_data <- 
  lapply(ensembl_datasets, function(ensembl_dataset){
    
    ## sometimes function useMart isn't able to connect to server and returns error
    # this chunck repeats useMart until there is no error
    mart <- "error"
    count <- 0
    while(class(mart) == "character"){
      count <- count + 1
      print(str_c(ensembl_dataset, " ", count))
      mart <- tryCatch(expr = useMart(biomart = ifelse(str_detect(ensembl_dataset, "mmusculus"), "ensembl", "ENSEMBL_MART_MOUSE"), dataset  = ensembl_dataset),
                       error = function(x) return("error"))
    }
    
    bm <- 
      getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "gene_exon_intron"),
            filters = c("ensembl_gene_id"),
            values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")],
            mart = mart) %>% 
      dplyr::mutate(ensembl_dataset = ensembl_dataset)
    
    return(bm)

  }) %>% 
  dplyr::bind_rows(.)

### write the sequences
if(!nrow(ensembl_data) == 0){
  
  # create dirs
  outdir <- file.path(outpath, gene_symbol)
  dir.create(outdir, showWarnings = F)
  
  # get and write sequences
  ensembl_seq <- 
    ensembl_data %>% 
    dplyr::filter(gene_exon_intron != "Sequence unavailable") %$% 
    gene_exon_intron %>% 
    DNAStringSet(.)
  names(ensembl_seq) <- 
    str_replace(ensembl_data$ensembl_dataset, "(_|\\.).+", "") %>% 
    str_c(., gene_symbol, sep = '_')
  
  ensembl_seq %<>% 
    .[order(width(.))] %>% 
    .[!duplicated(names(.))] %T>% 
    writeXStringSet(., file.path(outdir, str_c(gene_symbol, length(.), "fa", sep = ".")), format = 'fasta')
  
}


### test
# get ranges of genes
ensembl_ranges <- 
  ensembl_data %>% 
  as.tibble(.) %>% 
  dplyr::select(-gene_exon_intron) %>% 
  # dplyr::mutate(chromosome_name = stringr::str_c("chr", chromosome_name), 
  #               strand = ifelse(strand == 1, "+", "-")) %>% 
  dplyr::select(seqnames = chromosome_name, start = start_position, end = end_position, strand, gene_id = ensembl_gene_id, ensembl_dataset)

ensembl_data2 <- 
  lapply(1:nrow(ensembl_ranges), function(row_n){
    
    ensembl_dataset <- ensembl_ranges[row_n, ]
    
    ## sometimes function useMart isn't able to connect to server and returns error
    # this chunck repeats useMart until there is no error
    mart <- "error"
    count <- 0
    while(class(mart) == "character"){
      count <- count + 1
      print(str_c(ensembl_dataset$ensembl_dataset, " ", count))
      mart <- tryCatch(expr = useMart(biomart = ifelse(str_detect(ensembl_dataset$ensembl_dataset, "mmusculus"), "ensembl", "ENSEMBL_MART_MOUSE"), 
                                      dataset = ensembl_dataset$ensembl_dataset),
                       error = function(x) return("error"))
    }
    
    bm <- 
      getBM(attributes = c("ensembl_gene_id", "gene_exon_intron"),
            filters = c("chromosome_name", "start", "end", "strand"),
            values = list(ensembl_dataset$seqnames, ensembl_dataset$start, ensembl_dataset$end, ensembl_dataset$strand),
            mart = mart) 
    
    return(bm)
    
  }) %>% 
  dplyr::bind_rows(.)

