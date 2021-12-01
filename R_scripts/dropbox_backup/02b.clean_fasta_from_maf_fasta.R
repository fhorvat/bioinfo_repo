### INFO: 
### DATE: Mon Jul 08 12:52:07 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/UCSC_maf")

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

library(rtracklayer)
library(Biostrings)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# lnc1 multiz60way fasta path
lnc1_multiz_path <- file.path(inpath, "lnc1_3prime_end.UCSC.multiz60way.20190708.mafasta")

# lnc1 mouse strains fasta path
lnc1_strains_path <- file.path(inpath, "lnc1_3prime_end.UCSC.mouse_strains.20190707.mafasta")

# strains order path
strains_order_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/genomes", 
                                "genome_info.mouse_strains_order.txt")

######################################################## READ DATA
# read lnc1 multiz60way fasta
lnc1_multiz <- Biostrings::readDNAStringSet(lnc1_multiz_path)

# read lnc1 mouse strains fasta
lnc1_strains <- Biostrings::readDNAStringSet(lnc1_strains_path)

# read mouse strains order
strains_order <- readr::read_lines(strains_order_path)

######################################################## MAIN CODE
### clean fasta and coordinates
# lnc1 multiz60way
lnc1_multiz_tidy <-
  tibble(seq_name = names(lnc1_multiz), 
         seq_fasta = as.character(lnc1_multiz) %>% unname(.)) %>% 
  tidyr::separate(seq_name, c("species", "other"), sep = "\\.") %>%
  dplyr::group_by(species) %>% 
  dplyr::summarise(seq_fasta = str_c(seq_fasta, collapse = "")) %>% 
  dplyr::filter(species == "rn5")

# lnc1 mouse strains 
lnc1_strains_tidy <-
  tibble(seq_name = names(lnc1_strains), 
         seq_fasta = as.character(lnc1_strains) %>% unname(.)) %>% 
  tidyr::separate(seq_name, c("species", "other"), sep = "\\.") %>%
  dplyr::group_by(species) %>% 
  dplyr::summarise(seq_fasta = str_c(seq_fasta, collapse = "")) %>% 
  dplyr::filter(species != "rn6")

# join, remove gaps
lnc1_all <- 
  bind_rows(lnc1_multiz_tidy, lnc1_strains_tidy) %>% 
  dplyr::mutate(seq_fasta = str_remove_all(seq_fasta, "-")) %>% 
  dplyr::mutate(species = factor(species, levels = strains_order)) %>% 
  dplyr::arrange(species)

# create DNAStringSet, save as fasta
lnc1_fasta <- 
  lnc1_all$seq_fasta %>% 
  set_names(lnc1_all$species) %>% 
  Biostrings::DNAStringSet(x = .) %T>%
  Biostrings::writeXStringSet(x = ., filepath = file.path(outpath, "lnc1_3prime_end.UCSC.mouse_strains_rat.20190708.fasta"))

# do MSA using ClustalW, write alignment
lnc1_fasta_msa <- msa::msaClustalW(inputSeqs = lnc1_fasta)

# write
lnc1_fasta_msa@unmasked %T>% 
  Biostrings::writeXStringSet(., file = file.path(outpath, "lnc1_3prime_end.UCSC.mouse_strains_rat.20190708.ClustalW.msa.fasta"))

# write as .pdf
msaPrettyPrint(lnc1_fasta_msa, output = "pdf", paperWidth = 23, paperHeight = 5, 
               logoColors = "chemical",
               file = "lnc1_3prime_end.UCSC.mouse_strains_rat.20190708.ClustalW.msa.pdf", askForOverwrite = F)



