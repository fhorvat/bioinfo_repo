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
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS
# gets gene ID from different assemblies 
GetIDsMouseStrains <- function(id){
  
  ### sometimes function useMart isn't able to connect to server and returns error
  # this chunck repeats useMart until there is no error
  mart <- "error"
  count <- 0
  while(class(mart) == "character"){
    count <- count + 1
    print(str_c("m129s1svimj_gene_ensembl", " ", count))
    mart <- tryCatch(expr = useMart("ENSEMBL_MART_MOUSE", dataset = "m129s1svimj_gene_ensembl"), 
                     error = function(x) return("error"))
  }
  
  ### because it is possible to fetch only 6 attributes in one call, getBM runs in loop
  l.bm <- list()
  for(i in unique(set)){
    print(i)
    l.bm[[i]] <- unlist(getBM(attributes = c(attrs[i == set]), filters = c("ensembl_gene_id"), values = id, mart = mart))
  }
  bm <- unlist(l.bm)
  names(bm) <- str_replace(names(bm),'_.+','')
  return(bm)
  
}

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
              str_replace_all(., "Mouse | genes.*", ""))

######################################################## MAIN CODE
# create list of homolog ensembl genes
attrs <- 
  str_c(animals, "_homolog_ensembl_gene") %>% 
  .[. != "m129s1svimj_homolog_ensembl_gene"]
set <- rep(1:((length(attrs) %/% 3) + 1), each = 3)[1:length(attrs)]
ensembl_datasets <- str_c(animals, "_gene_ensembl", sep = "")

### gets sequence of gene from all animals
# create dirs
gene_symbol <- "Ago2"
gene_id <- "MGP_129S1SvImJ_G0021919"
outdir <- file.path(outpath, gene_symbol)
dir.create(outdir, showWarnings = F)

# get geneIDs
ids <- 
  GetIDsMouseStrains(id = gene_id) %>% 
  .[!duplicated(.)]

### loops over the organisms and fetches sequences
d.seq <- list()
for(ensembl_dataset in ensembl_datasets){
  
  # sometimes function useMart isn't able to connect to server and returns error
  # this chunck repeats useMart until there is no error
  mart <- "error"
  count <- 0
  while(class(mart) == "character"){
    count <- count + 1
    print(str_c(ensembl_dataset, " ", count))
    mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_MOUSE", dataset  = ensembl_dataset),
                     error = function(x) return("error"))
  }
  
  bm <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "gene_exon_intron"),
              filters = c("ensembl_gene_id"),
              values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")],
              mart = mart)
  d.seq[[ensembl_dataset]] <- bm
  
}

### writes the sequences
if(!nrow(d.seq) == 0){
  
  d.seq <- d.seq[d.seq$sequence != "Sequence unavailable", ]
  s <- DNAStringSet(d.seq$sequence)
  names(s) <- str_replace(rownames(d.seq), "(_|\\.).+", "")
  s <- s[order(-width(s))]
  s <- s[!duplicated(names(s))]
  names(s) <- str_c(names(s), gene_id, sep = '_')
  
  writeXStringSet(s, file.path(outdir, str_c(gene_id, length(s), "fa", sep = ".")), format = 'fasta')
}



