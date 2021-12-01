### INFO: 
### DATE: Fri Jul 05 22:51:19 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/Muridae")

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

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# other rodents genomes info path
rodents_genomes_info_path <- file.path(inpath, "genomes_info.Muridae.csv")

######################################################## READ DATA
# read other rodents genomes info 
rodents_genomes_info <- readr::read_csv(rodents_genomes_info_path)

######################################################## MAIN CODE
# clean other rodents genomes info
rodents_genomes_info_tidy <- 
  rodents_genomes_info %>% 
  dplyr::mutate(assembly_name = str_remove(GenBank_accession, "GCA_.* "), 
                GenBank_accession = str_remove(GenBank_accession, assembly_name) %>% str_trim(.), 
                GenBank_URL = str_c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/", 
                                    str_replace_all((GenBank_accession %>% str_remove_all(., "^GCA_|\\.[1-9]{1,}$")), "(.{3})", "\\1/"),
                                    GenBank_accession, "_", assembly_name, "/",
                                    GenBank_accession, "_", assembly_name, "_genomic.fna.gz")) %>% 
  dplyr::mutate(genome_name_short = str_remove_all(organism_scientific, "(?<=.{3}).*(?= )|(?<= .{3}).*") %>% str_to_title(.) %>% str_remove_all(., " "), 
                genome_name_short = str_replace_all(genome_name_short, c("MusMus" = "mm10", "RatNor" = "rn6"))) %>% 
  dplyr::select(organism_scientific, assembly_name, genome_name_short, organism_common, GenBank_accession, UCSC_URL, GenBank_URL) %>% 
  arrange(organism_scientific) %T>% 
  readr::write_csv(., file.path(outpath, "genomes_info.Muridae.tidy.csv"))
