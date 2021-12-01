#!/home/students/fhorvat/R/bin/Rscript
### INFO: reads sequencing table, writes sample table and links to raw files
### DATE: 16. 10. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/documentation/sample_table")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(readxl)

######################################################## PATH VARIABLES
outpath <- getwd()
raw_path <- "/common/RAW/Svoboda/Split_RNAseq_2017"
links_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
original_table <- readxl::read_excel(path = list.files(outpath, pattern = "*.xlsx", full.names = T))

######################################################## MAIN CODE
# list of fastq files
sample_table <- 
  tibble(file_path = list.files(path = raw_path, pattern = "*txt.gz", full.names = T, recursive = T)) %>% 
  dplyr::mutate(file_name = (file_path %>% basename()),
                sample_name = (stringr::str_extract(file_name, "BC[0-9]{1,2}.*_|SG.*_") %>% stringr::str_replace(., "_17.*", "")), 
                seq_name = (stringr::str_extract(file_name, "17.*") %>% stringr::str_replace(., "-1-1_Malik.*", ""))) %>% 
  dplyr::select(-file_name)

# output table
out_table <- 
  original_table %>% 
  dplyr::mutate(sample_name = stringr::str_replace_all(`Sample NGS ID`, " ", "_")) %>% 
  dplyr::left_join(sample_table, by = "sample_name") %>% 
  dplyr::mutate(organism = str_replace(Sample, " [[:digit:]]", "") %>% tolower(.)) %>% 
  dplyr::mutate(organism = replace(organism, organism %in% c("cast", "pwd"), "mouse_strains"), 
                organism = replace(organism, !(organism %in% c("cow", "rat", "rabbit", "mouse_strains")), "mouse"), 
                resequnced = ifelse(str_detect(file_path, "CBH5JACXX"), ".RS", ""), 
                sample_name = str_c(sample_name, resequnced)) %>%
  dplyr::select(-resequnced) %>% 
  dplyr::select(`Sample NGS ID`, sample_name,  seq_name, organism, everything()) %T>%
  readr::write_csv(., path = file.path(outpath, "sample_table.csv"))

# links to single-end fastq files
make_links_single_end <- 
  out_table %>% 
  dplyr::filter(str_detect(string = `sequencing type`, pattern = "SE")) %>% 
  dplyr::mutate(link_name = str_c("s_", sample_name, ".SE.txt.gz"), 
                make_links = (str_c("ln -s ", file_path, " ",
                                    file.path(links_path, "mouse/Data/Raw/Links", link_name)))) %$%
  make_links %T>% 
  readr::write_lines(., file.path(outpath, "make_links.SE.sh"))

# links to pair-end fastq files
make_links_pair_end <-
  out_table %>%
  dplyr::filter(str_detect(string = `sequencing type`, pattern = "PE")) %>%
  dplyr::mutate(link_name = str_c("s_", str_replace(sample_name, "Cow[0-9]{1}_", ""), "_",
                                  str_replace(Sample, " ", ""), ".PE_",
                                  rep(1:2, times = (length(file_path) / 2)), ".txt.gz"),
                make_links = (str_c("ln -s ", file_path, " ",
                                    file.path(links_path, "paired_end", organism, "Data/Raw/Links", link_name)))) %$%
  make_links %T>%
  readr::write_lines(., file.path(outpath, "make_links.PE.sh"))
