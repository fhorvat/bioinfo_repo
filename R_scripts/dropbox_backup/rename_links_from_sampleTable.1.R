### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Raw/Links")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# inpath
inpath <- getwd()

# sample table
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/documentation/RNAseq_2016_11_23_sampleTable.csv"

# raw files
raw_path <- list.files(path = "/common/RAW/Svoboda/2016-11-23-HT3LFBGXY", pattern = "*.txt.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "SE"

# create rename executable
rename_file <- 
  sample_table %>% 
  dplyr::select(sample_id = sample, genotype = Type, stage = Sample, barcode = `most common barcode`) %>% 
  dplyr::mutate(genotype = stringr::str_remove(genotype, " B6"), 
                name = str_c("s", genotype, sep = "_"), 
                resequencing = str_detect(sample_id, "HVF5KBGXY")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(barcode_n = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sequence = ifelse(), 
                name = str_c(name, "_r", sequence),
                name = ifelse(resequencing, str_c(name, ".reseq"), name), 
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>% 
  dplyr::select(sample_id = name, stage, genotype) %T>% 
  readr:::write_csv(., path = file.path(sra_path, str_c(basename(sra_path), ".sampleTable.csv")))

# only for PE data
if(pairing == "PE"){
  
  rename_file %<>% 
    rbind(., .) %>% 
    dplyr::group_by(run) %>% 
    dplyr::mutate(sequence = 1:n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(run = str_c(run, "_", sequence), 
                  name = str_c(name, "_", sequence)) %>% 
    dplyr::select(-sequence)
  
}


# save make links script
rename_file %<>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = "make_links.sh")

