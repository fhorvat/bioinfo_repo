### INFO: 
### DATE: Wed Oct 03 16:20:53 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/ensembl_counts")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)

library(biomaRt)
library(googlesheets)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome info path
genome_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genome_info.csv"

######################################################## READ DATA
# # load google sheets token and read list of genomes from GSheets, write genomes to lobsang
# suppressMessages(gs_auth(token = "/common/WORK/fhorvat/code_library/R_scripts/googlesheets_token.rds", verbose = FALSE))
# genome_df <- gs_read(gs_title("genome_assemblies"))
# readr::write_csv(genome_df, "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genome_info.csv")

# read genomes info
genome_df <- readr::read_csv(genome_info_path)

######################################################## MAIN CODE
# get ensembl names - separate mouse assemblies from other
animal_ensembl <- 
  genome_df %>% 
  dplyr::filter(!(name %in% c("mouse_cast", "mouse_PWD"))) %$% 
  ensembl_name

# mouse names
mouse_ensembl <- 
  genome_df %>% 
  dplyr::filter(name %in% c("mouse_cast", "mouse_PWD")) %$% 
  ensembl_name

# set ensembl name and version
ensembl_name <- "mmusculus"
ensembl_release <- 93

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(93, 92, 91, 89, 86),
         date = c("Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_release) %$%
  URL_archive

### load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  
  count <- count + 1
  print(str_c(ensembl_name, "_gene_ensembl", " ", count))
  
  # load ENSEMBL mart
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
                   error = function(x) return("error"))
  
  # stop if count get too big
  if(count > 2){
    stop("Something's not right")
  }
  
}

# load homologs od mmusculus genes from animal datasets - in loop
animals_homologs <- 
  purrr::map(animal_ensembl[animal_ensembl != "mmusculus"], function(animal){
    
    # load mmusculus homolog ensembl gene
    getBM(attributes = c("ensembl_gene_id", str_c(animal, "_homolog_ensembl_gene")), mart = mart) %>% 
      as.tibble(.) %>% 
      dplyr::filter((!!sym(str_c(animal, "_homolog_ensembl_gene"))) != "") %>% 
      group_by(!!sym(str_c(animal, "_homolog_ensembl_gene"))) %>% 
      filter(n() == 1) %>% 
      dplyr::ungroup()


  })

# get homologs of mmusculus genes in mouse strains
mouse_strains_homologs <- 
  purrr::map(mouse_ensembl, function(mouse_strain){
    
    # load mmusculus homolog ensembl gene
    getBM(attributes = c("mmusculus_homolog_ensembl_gene", "ensembl_gene_id"), 
               mart = useMart(biomart = "ENSEMBL_MART_MOUSE", dataset = str_c(mouse_strain, "_gene_ensembl"))) %>% 
      as.tibble(.) %>% 
      dplyr::filter(mmusculus_homolog_ensembl_gene != "") %>% 
      magrittr::set_colnames(., c("ensembl_gene_id", str_c(mouse_strain, "_homolog_ensembl_gene")))
    
  })

# join
all_homologs <- c(animals_homologs, mouse_strains_homologs)

# reduce join
x <- purrr::reduce(all_homologs, left_join, by = "ensembl_gene_id")


