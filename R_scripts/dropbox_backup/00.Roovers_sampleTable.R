### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Roovers_2015_CellRep_GSE64942"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq", reference, "Data/Raw/Links"))

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# SRA path
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/smallRNAseq", reference)

# rename file path
rename_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Roovers_2015_CellRep_GSE64942/Data/Raw/Links/02_rename_samples.sh"

######################################################## READ DATA
# read rename file
sample_tb <- readr::read_delim(rename_path, delim = " ", col_names = c("command", "run", "sample_id"))

######################################################## MAIN CODE
# clean sample table
sample_tb_clean <- 
  sample_tb %>% 
  dplyr::select(sample_id, run) %>% 
  dplyr::mutate(run = run %>% str_remove(., "\\.fastq\\.gz"),
                sample_id = sample_id %>% str_remove(., "\\.fastq\\.gz") %>% str_c(., ".SE"), 
                animal = sample_id %>% str_remove(., "^s_") %>% str_remove(., "_.*"), 
                oxidation = sample_id %>% str_extract(., "NaIO4") %>% str_replace_na(.) %>% str_replace(., "NA", "none"), 
                stage = sample_id %>% str_remove(., "^s_") %>% str_remove(., animal) %>% str_remove(., oxidation) %>% 
                  str_remove_all(., "_r[0-9]+|\\.SE$") %>% str_remove_all(., "^_|_$")) %>% 
  dplyr::select(sample_id, stage, animal, oxidation, run)

# save sample table
sample_tb_clean %T>%
  readr:::write_csv(., path = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., path = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

# save make links script
sample_tb_clean %>% 
  dplyr::arrange(sample_id) %>% 
  dplyr::mutate(name = str_c(outpath, "/", sample_id, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
