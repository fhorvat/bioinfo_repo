### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133")

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  
  count <- count + 1
  print("mmusculus_gene_ensembl", " ", count)
  
  # load ENSEMBL mart
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "http://jul2018.archive.ensembl.org"),
                   error = function(x) return("error"))
  
  # stop if count get too big
  if(count > 2){
    stop("Something's not right")
  }
  
}

# get info about genes
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type"), mart = mart) %>% 
  as.tibble(.) %>% 
  dplyr::filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::rename(gene_id = ensembl_gene_id) %T>% 
  readr::write_csv(., file.path(outpath, "ensembl.93.mouse_human.one2one_homologs.csv"))
