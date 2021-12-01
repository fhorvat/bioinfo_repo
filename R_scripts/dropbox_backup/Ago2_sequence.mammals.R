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
# gets gene ID from different animals/assemblies 
GetIDs <- function(id, orgs){
  
  # create list of homolog ensembl genes
  orgs <- orgs[names(orgs) != "human"]
  attrs <- stringr::str_c(orgs, "_homolog_ensembl_gene", sep = "")

  # sometimes function useMart isn't able to connect to server and returns error
  # this chunck repeats useMart until there is no error
  mart <- "error"
  count <- 0
  while(class(mart) == "character"){
    count <- count + 1
    print(str_c("hsapiens_gene_ensembl", " ", count))
    mart <- tryCatch(expr = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                     error = function(x) return("error"))
  }

  # get Ensembl gene ID from gene symbol
  hbm <- getBM(attributes = c("ensembl_gene_id"), filters = c("hgnc_symbol"), values = id, mart = mart)
  names(hbm) <- "hsapiens"
  
  ### because it is possible to fetch only 6 attributes in one call, getBM runs in loop
  # get length of attributes	
  l <- length(attrs)
  set <- rep(1:((l %/% 3) + 1), each = 3)[1:l]
  l.bm <- list()

  # get Ensembl gene IDs of homolog gene for all other animals/assemblies
  for(i in unique(set)){
    print(i)
    l.bm[[i]] <- unlist(getBM(attributes = attrs, filters = c("hgnc_symbol"), values = c(id), mart = mart))
  }
  bm <- unlist(l.bm)
  names(bm) <- str_replace(names(bm),'_.+','')
  return(bm)
}

######################################################## READ DATA

######################################################## MAIN CODE
# set animals from which you want sequences 
animals <- c('guinea pig' = 'cporcellus',
             'rat' = 'rnorvegicus',
             'mouse' = 'mmusculus',
             'rabbit' = 'ocuniculus',
             # 'squirrle' = 'itridecemlineatus',
             'kangaroo rat' = 'dordii',
             'human' = 'hsapiens',
             'pig' = 'sscrofa'
             # 'elephant' = 'lafricana',
             # 'bat' = 'mlucifugus',
             # 'platypus' = 'oanatinus',
             # 'dog' = 'cfamiliaris',
             # 'cow' = 'btaurus'
)

# check availible ensembl datasets, filter animals
available_ensembl_datasets <- 
  useMart(biomart = "ensembl") %>% 
  listDatasets(.) %$%
  dataset
ensembl_datasets <- str_c(animals, "_gene_ensembl", sep = "")
ensembl_datasets <- ensembl_datasets[ensembl_datasets %in% available_ensembl_datasets]
  
# set genes you want sequences
hgnc <- c("Ago2") 

### gets sequence of gene from all animals
# create dirs
outdir <- file.path(outpath, hgnc)
dir.create(outdir, showWarnings = F)

# get geneIDs
ids <- GetIDs(hgnc, animals)
ids <- ids[!duplicated(ids)]

### loops over the organisms and fetches sequences
l.data <- list()
for(ensembl_dataset in ensembl_datasets){
  
  # sometimes function useMart isn't able to connect to server and returns error
  # this chunck repeats useMart until there is no error
  mart <- "error"
  count <- 0
  while(class(mart) == "character"){
    count <- count + 1
    print(str_c(ensembl_dataset, " ", count))
    mart <- tryCatch(expr = useMart(biomart = "ensembl", dataset  = ensembl_dataset),
                     error = function(x) return("error"))
  }
  
  sequence <- getSequence(id = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")], 
                          type = "ensembl_gene_id",
                          seqType = "gene_exon_intron", 
                          mart = mart)
  
  # bm <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "coding"), 
  #             filters = c("ensembl_gene_id"),
  #             values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")], 
  #             mart = mart)
  l.data[[ensembl_dataset]] <- sequence

}

# bind to data frame
d.seq <- 
  do.call(rbind, l.data) %>% 
  data.table::setnames(., old = 1, new = "sequence")

### writes the sequences
if(!nrow(d.seq) == 0){
  
  d.seq <- d.seq[d.seq$sequence != "Sequence unavailable", ]
  s <- DNAStringSet(d.seq$sequence)
  names(s) <- str_replace(rownames(d.seq), "(_|\\.).+", "")
  s <- s[order(-width(s))]
  s <- s[!duplicated(names(s))]
  names(s) <- str_c(names(s), hgnc, sep = '_')
  
  writeXStringSet(s, file.path(outdir, str_c(hgnc, length(s), "fa", sep = ".")), format = 'fasta')
}



