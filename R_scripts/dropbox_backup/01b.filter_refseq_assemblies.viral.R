### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/other/bacterial_contamination/genomes/viral")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# splits string into equal-length parts
str_dice <- function(s, width) { 
  L <- str_length(s)
  str_sub(s, start = seq(1L, L - width + 1, width), end = seq(width, L, width))
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get list of genomes
genomes_path <- file.path(inpath, "genomes.txt")

######################################################## READ DATA
# read genomes list
genomes_tb <- readr::read_delim(genomes_path, delim = "\t", col_names = c("assembly_name", "full_species", "refseq_accession"))

######################################################## MAIN CODE
# tidy dataframe, get one assembly per species, construct FTP links
genomes_tidy <- 
  genomes_tb %>% 
  dplyr::mutate(genus = str_remove(full_species, " .*"), 
                species = str_remove(full_species, "\\w+ ") %>% str_remove(., " .*"), 
                genus_species = str_c(genus, " ", species)) %>% 
  dplyr::group_by(genus_species) %>% 
  dplyr::slice_sample(n = 1)

# create links
genome_links <- 
  genomes_tidy %>% 
  ungroup(.) %>% 
  dplyr::select(refseq_accession, assembly_name) %>% 
  dplyr::rowwise(.) %>%  
  dplyr::mutate(accession_link = str_remove(refseq_accession, "GCF_") %>% 
                  str_remove(., "\\..*") %>% 
                  str_dice(., 3) %>% 
                  str_c(., collapse = "/")) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(full_link = str_c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/", 
                                  accession_link, "/",
                                  refseq_accession, "_", assembly_name, "/", 
                                  refseq_accession, "_", assembly_name, "_genomic.fna.gz")) %$%
  full_link %T>% 
  readr::write_lines(., file.path(outpath, "links.txt"))
