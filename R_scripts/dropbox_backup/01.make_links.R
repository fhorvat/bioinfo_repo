### INFO: 
### DATE: Fri Jul 05 22:51:19 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/genomes")

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

# mouse strains path
mouse_genomes_path <- "/common/DB/genome_reference/mouse/strains"

# mouse strains genomes info path
mouse_genomes_info_path <- file.path(mouse_genomes_path, "genomes_info.mouse_strains.tidy.csv")

# other rodents path
rodents_genomes_path <- "/common/DB/genome_reference/Muridae"

# other rodents genomes info path
rodents_genomes_info_path <- file.path(rodents_genomes_path, "genomes_info.Muridae.tidy.csv")

######################################################## READ DATA
# read mouse genomes info 
mouse_genomes_info <- readr::read_csv(mouse_genomes_info_path)

# read other rodents genomes info 
rodents_genomes_info <- readr::read_csv(rodents_genomes_info_path)

######################################################## MAIN CODE
# join and save
rodents_genomes_info <- 
  bind_rows(mouse_genomes_info, rodents_genomes_info) %>% 
  arrange(organism_scientific) %>% 
  dplyr::mutate(URL = ifelse(!is.na(UCSC_URL), UCSC_URL, GenBank_URL)) %>% 
  dplyr::select(-UCSC_URL, -GenBank_URL) %T>% 
  readr::write_csv(., file.path(outpath, "genome_info.Muridae_with_mouse_strains.tidy.20190706.csv"))

# rename links
links_fasta <- 
  tibble(genome_file = list.files(c(mouse_genomes_path, rodents_genomes_path), pattern = "\\.fa\\.gz$", full.names = T), 
         genome_name = genome_file %>% basename(.)) %>%
  left_join(., rodents_genomes_info %>% 
              dplyr::mutate(genome_name = str_c(GenBank_accession, "_", assembly_name, ".fa.gz")) %>% 
              dplyr::select(genome_name, genome_name_short), 
            by = "genome_name") %>% 
  dplyr::mutate(genome_new_file = str_c(genome_name_short, ".fa.gz"), 
                link = str_c("ln -s", 
                             genome_file,
                             file.path(outpath, genome_new_file), sep = " ")) %$%
  link 

# add index files to links
links_fai <- links_fasta %>% str_replace_all(., ".fa.gz", ".fa.gz.fai")
links_gzi <- links_fasta %>% str_replace_all(., "\\.fa\\.gz", ".fa.gz.gzi")

# write renaming script
readr::write_lines(c(links_fasta, links_fai, links_gzi), file.path(outpath, "make_links.sh"))
