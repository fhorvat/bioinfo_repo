### INFO: 
### DATE: Fri Jul 05 22:51:19 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/mouse/strains")

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

# mouse strains genomes info path
mouse_genomes_info_path <- file.path(inpath, "genomes_info.mouse_strains.txt")

######################################################## READ DATA
# read mouse genomes info 
mouse_genomes_info <- readr::read_lines(mouse_genomes_info_path)

######################################################## MAIN CODE
# clean mouse strains genome info
mouse_genomes_info_tidy <- 
  tibble(raw = mouse_genomes_info) %>% 
  dplyr::filter(raw != "") %>% 
  tidyr::separate(raw, into = c("category", "value"), sep = "\\s", extra = "merge") %>% 
  dplyr::filter(category %in% c("genome", "twoBitPath", "organism")) %>% 
  dplyr::mutate(species = rep(1:(nrow(.) / 3), each = 3)) %>% 
  tidyr::spread(category, value) %>% 
  dplyr::select(-species) %>% 
  dplyr::mutate(assembly_name = str_c(genome, "_v1"), 
                organism_common = ifelse(organism == "house mouse", str_c(organism, assembly_name %>% str_replace(., "_", "/"), sep = " "), organism), 
                GenBank_accession = twoBitPath %>% str_remove(., str_c("_", assembly_name, ".*")), 
                UCSC_URL = file.path("http://hgdownload.soe.ucsc.edu/hubs/mouseStrains", twoBitPath), 
                GenBank_URL = str_c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/", 
                                    str_replace_all((GenBank_accession %>% str_remove_all(., "^GCA_|\\.1$")), "(.{3})", "\\1/"),
                                    str_replace(twoBitPath, "\\.2bit$", "_genomic.fna.gz"))) %>%
  dplyr::mutate(organism_scientific = "Mus musculus",
                organism_scientific = replace(organism_scientific, assembly_name == "CAST_EiJ_v1", "Mus musculus castaneus"),
                organism_scientific = replace(organism_scientific, assembly_name == "PWK_PhJ_v1", "Mus musculus musculus"),
                organism_scientific = replace(organism_scientific, assembly_name == "WSB_EiJ_v1", "Mus musculus domesticus"), 
                organism_scientific = replace(organism_scientific, assembly_name == "SPRET_EiJ_v1", "Mus spretus")) %>% 
  dplyr::mutate(genome_name_short = assembly_name) %>% 
  dplyr::select(organism_scientific, assembly_name, genome_name_short, organism_common, GenBank_accession, UCSC_URL, GenBank_URL) %>% 
  arrange(organism_scientific) %T>% 
  readr::write_csv(., file.path(outpath, "genomes_info.mouse_strains.tidy.csv"))
