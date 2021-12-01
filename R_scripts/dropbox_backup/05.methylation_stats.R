### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Analysis")

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

library(Biostrings)
library(BSgenome.Maur.UCSC.Siomi)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# genome fasta file path
fasta_path <- file.path(genome_path, "hamster.sequel.draft-20200302.arrow.fasta")

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi.trimmed"

# coverage files path
coverage_path <- list.files(mapped_path, pattern = ".*deduplicated\\..*\\.bw", full.names = T)

# raw coverage path
raw_coverage_path <- coverage_path[str_detect(basename(coverage_path), "\\.raw\\.bw")]

# metyhlation coverage path
meth_coverage_path <- coverage_path[!(str_detect(basename(coverage_path), "\\.raw\\.bw"))]

######################################################## READ DATA

######################################################## MAIN CODE
### total genome coverage
# read raw coverage data, calculate total genome coverage
coverage_tb_full <- purrr::map(raw_coverage_path, function(path){
  
  # read bw
  bw_gr <- rtracklayer::import.bw(path, as = "RleList")
  
  # get lengths of scaffolds
  scaffolds_tb <- 
    lapply(bw_gr, length) %>% 
    unlist(.)
  
  # create table
  scaffolds_tb <- tibble(scaffold = names(scaffolds_tb), 
                         scaffold_length = unname(scaffolds_tb))
  
  # get coverage length
  scaffold_coverage_tb <- 
    lapply(bw_gr, function(x) x[x > 0] %>% length(.)) %>%  
    unlist(.)
  
  # create table
  scaffolds_coverage_tb <- tibble(scaffold = names(scaffold_coverage_tb), 
                                  scaffold_coverage = unname(scaffold_coverage_tb))
  
  # join together, calculate percentage
  coverage_tb <-
    left_join(scaffolds_tb, scaffolds_coverage_tb, by = "scaffold") %>% 
    dplyr::summarise(scaffold_length = sum(scaffold_length), 
                     scaffold_coverage = sum(scaffold_coverage)) %>% 
    dplyr::mutate(coverage_perc = 100*(scaffold_coverage / scaffold_length),
                  coverage_perc = round(coverage_perc, 2), 
                  sample_id = basename(path) %>% str_remove("_bismark_bt2_pe.deduplicated.raw.bw")) %>% 
    dplyr::select(sample_id, raw_coverage_percentage = coverage_perc)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)


### methylation coverage
# read CpG methylation coverage
methylation_coverage_tb <- purrr::map(meth_coverage_path, function(path){
  
  # read bw
  bw_gr <- rtracklayer::import.bw(path, as = "RleList")
  
  # get coverage length
  cpg_tb <- 
    lapply(bw_gr, function(x) x[x > 0] %>% length(.)) %>%  
    unlist(.)
  
  # create table
  cpg_tb <- 
    tibble(scaffold = names(cpg_tb), 
           cpg_number = unname(cpg_tb)) %>% 
    dplyr::summarise(cpg_number = sum(cpg_number)) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove("_bismark_bt2_pe.deduplicated.bw")) %>% 
    dplyr::select(sample_id, cpg_number)
  
  # return
  return(cpg_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# # get the total number of CpG in fasta file
# cpg_total <- system(command = str_c('grep -v "^>" ', fasta_path, '  | tr -d "\n" | tr -c "[GCgc]" "\n" | grep -v "^$" | grep -o -i "CG" | wc -l'),
#                     intern = TRUE) 
# cpg_total <- as.numeric(cpg_total)

# get the total number of CpG in BSGenome
cpg_total <- (length(vmatchPattern("CG", BSgenome.Maur.UCSC.Siomi)) / 2)

# add the total number of CpGs to the table
methylation_tb <- 
  methylation_coverage_tb %>% 
  dplyr::mutate(cpg_total = cpg_total,
                cpg_perc = 100*(cpg_number / cpg_total),
                cpg_perc = round(cpg_perc, 2)) %>% 
  dplyr::select(sample_id, cpg_perc)

# join to one table, save
out_tb <- 
  coverage_tb_full %>% 
  dplyr::left_join(., methylation_tb, by = "sample_id") %T>%
  readr::write_csv(., file.path(outpath, str_c("hamster_oocyte_Mov10l.bisulfite", 
                                               "Bismark_Siomi.trimmed.methylations_stats.csv", 
                                               sep = ".")))  





