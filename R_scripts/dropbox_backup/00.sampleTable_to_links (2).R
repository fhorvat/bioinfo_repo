### INFO: 
### DATE: Fri Sep 28 08:18:54 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Raw/Links")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# documentation path 
documentation_path <- file.path(getwd(), "../../Documentation")

# sample table
sample_table_path <- list.files(path = documentation_path, ".*sampleTable.raw.csv", full.names = T)

# raw files
raw_path <- list.files(path = "/common/RAW/Svoboda/smallRNASeq_2018/2018-09-28-CCME6ANXX", pattern = "*.txt.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "SE"

# raw files table
raw_table <- 
  tibble(raw_path) %>% 
  dplyr::mutate(sample_n = raw_path %>% 
                  basename(.) %>% 
                  str_remove_all(., "CCME6ANXX_TD-dsRNA_18s004566-1-1_Malik_lane4|_sequence.txt.gz"),
                replicate = stringr::str_extract(sample_n, "II$|I$"),
                replicate = ifelse(replicate == "II", 2, 1), 
                sample_n = str_remove(sample_n, "II$|I$"))
                
# create rename executable
rename_file <- 
  sample_table %>% 
  dplyr::select(sample_n = `Sample No`, transfection = Transfection, cell = `cell type`) %>% 
  dplyr::mutate(sample_n = as.character(sample_n), 
                transfection = stringr::str_replace(transfection, "-", "_"), 
                plasmid = stringr::str_extract(transfection, "pCag|pU6"),
                genotype = stringr::str_extract(transfection, "Mos|Lin28|Elav2"), 
                plasmid = ifelse(str_detect(transfection, "long"), str_c(plasmid, "_long"), plasmid), 
                name = str_c("s_", transfection)) %>% 
  dplyr::bind_rows(., .) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_c(name, "_r", replicate, ".", pairing)) %>%
  dplyr::left_join(., raw_table, by = c("sample_n", "replicate")) %>% 
  dplyr::select(sample_id, sample_n, transfection, plasmid, genotype, cell, replicate, raw_path) %>% 
  dplyr::arrange(sample_n)

# save sample table
rename_file %T>% 
  readr:::write_csv(., path = file.path(documentation_path, str_replace(basename(sample_table_path), "sampleTable.raw.csv", "sampleTable.csv")))

rename_file %T>% 
  readr:::write_csv(., path = file.path(unique(dirname(raw_path)), str_replace(basename(sample_table_path), "sampleTable.raw.csv", "sampleTable.csv")))
      
# save make links script
rename_file %>% 
  dplyr::mutate(name = str_c(outpath, "/", sample_id, ".txt.gz"), 
                make_links = str_c("ln -s ", raw_path, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = "make_links.sh")
