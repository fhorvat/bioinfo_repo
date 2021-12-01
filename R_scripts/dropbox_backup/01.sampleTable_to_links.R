### INFO: create renamed links from sample table
### DATE: Sat Jul 07 14:25:24 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Zuzka_3T3_PAPD7/Data/Raw/Links")

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

# documentation path 
documentation_path <- file.path(getwd(), "../../Documentation")

# sample table
sample_table_path <- list.files(path = documentation_path, ".*sampleTable.raw.csv", full.names = T)

# raw files
raw_path <- list.files(path = "/common/RAW/Svoboda/smallRNASeq_2018/2018-07-05-CCGH0ANXX/Zuzka_3T3_PAPD7", pattern = "*.txt.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "SE"

# raw files table
raw_table <- 
  tibble(raw_path) %>% 
  dplyr::mutate(sample_id = raw_path %>% 
                  basename(.) %>% 
                  str_remove_all(., "CCGH0ANXX_P7_3T3_3D6_18s003220-1-1_Malik_lane6|_sequence.txt.gz"))

# create rename executable
rename_file <- 
  sample_table %>% 
  dplyr::select(sample_id = `Sample name`, genotype = Genotype, transfection = Transfection, barcode = BARCODE) %>% 
  dplyr::mutate(genotype = stringr::str_replace(genotype, " ", "_"), 
                name = str_c("s", sample_id, sep = "_") %>% str_remove(., "[0-9]$"), 
                transfection = str_replace(transfection, "no transfection", "no_transfection") %>% 
                  str_replace(., " variant of ", "_")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(name = str_c(name, "_r", sequence, ".", pairing)) %>%
  dplyr::select(-sequence) %>% 
  dplyr::left_join(., raw_table, by = "sample_id")

# save sample table
rename_file %>% 
  dplyr::select(sample_id = name, sample = sample_id, genotype, transfection, raw_path) %T>% 
  readr:::write_csv(., path = file.path(documentation_path, str_replace(basename(sample_table_path), "sampleTable.raw.csv", "sampleTable.csv")))

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
rename_file %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"), 
                make_links = str_c("ln -s ", raw_path, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = "make_links.sh")

