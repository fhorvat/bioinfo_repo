### INFO: 
### DATE: Tue Jul 09 19:28:50 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/results")

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

# list fasta files
lnc1_fasta_paths <- list.files(inpath, ".*lnc1_3prime\\.fa", full.names = T)

# genomes path
genomes_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/genomes"

# species order path
species_order_path <- file.path(genomes_path, "genome_info.Muridae_order.txt")

# genomes info path
genomes_info_path <- file.path(inpath, "..", "lnc1_3prime.Muridae.genomes_info.csv")

######################################################## READ DATA
# read fasta
lnc1_seq <- 
  Biostrings::readDNAStringSet(lnc1_fasta_paths) %>% 
  reverseComplement(.)

# read species order
species_order <- readr::read_lines(species_order_path)

# read genomes info
genomes_info <- readr::read_csv(genomes_info_path)

######################################################## MAIN CODE
# add species to names
names(lnc1_seq) <- str_c(basename(lnc1_fasta_paths) %>% str_remove("\\.fa"), 
                         names(lnc1_seq), sep = ".")

# order and write as one file fasta
lnc1_seq_ordered <- 
  lnc1_seq %>% 
  .[names(.) %>% str_remove(., "\\.lnc1_3prime.*") %>% match(., species_order) %>% order(.)] %T>% 
  Biostrings::writeXStringSet(., "lnc1_3prime.Muridae.fasta")

# write as individual fasta
purrr::map(1:length(lnc1_seq_ordered), function(n){
  
  # get one sequence
  lnc1_seq_single <- lnc1_seq_ordered[n]
  
  # get name
  lnc1_seq_single_name <- 
    names(lnc1_seq_single) %>% 
    str_remove(., "\\.lnc1_3prime.*")
    
  # save as single .fasta
  Biostrings::writeXStringSet(lnc1_seq_single, str_c(lnc1_seq_single_name, ".lnc1_3prime.stranded.fa"))
  
})


### MSA
# do MSA using ClustalW
lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_seq)

# write ordered
lnc1_seq_ordered_msa <- 
  lnc1_msa@unmasked %>% 
  .[names(.) %>% str_remove(., "\\.lnc1_3prime.*") %>% match(., species_order) %>% order(.)] %T>% 
  Biostrings::writeXStringSet(., file = file.path(outpath, "lnc1_3prime.Muridae.ClustalW.msa.fasta"))

# write as .pdf
msaPrettyPrint(lnc1_msa, subset = lnc1_msa@unmasked %>% names(.) %>% str_remove(., "\\.lnc1_3prime.*") %>% match(., species_order) %>% order(.), 
               output = "pdf", paperWidth = 32, paperHeight = 5, 
               logoColors = "chemical",
               file = "lnc1_3prime.Muridae.ClustalW.msa.pdf", askForOverwrite = F)


### add to genomes table
# sequences to table
lnc1_seq_tb <- 
  tibble(species = names(lnc1_seq_ordered), 
         sequence = as.character(lnc1_seq_ordered)) %>% 
  dplyr::mutate(species = str_remove(species, "\\.lnc1_3prime.*"))

# MSA sequences to table
lnc1_seq_msa_tb <- 
  tibble(species = names(lnc1_seq_ordered_msa), 
         sequence_ClustalW = as.character(lnc1_seq_ordered_msa)) %>% 
  dplyr::mutate(species = str_remove(species, "\\.lnc1_3prime.*"))

# add all to one table, join with genomes info
lnc1_seq_tb %>% 
  dplyr::left_join(., genomes_info, by = "species") %>% 
  dplyr::left_join(., lnc1_seq_msa_tb, by = "species") %>% 
  dplyr::select(species, coordinates:organism_common, sequence, sequence_ClustalW) %T>%
  readr::write_csv(., "lnc1_3prime.Muridae.genomes_info.with_sequences.csv")

