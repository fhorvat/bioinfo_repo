### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/2019_Dec/Data/Documentation")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw files path
raw_path <- "/common/RAW/Svoboda/Papd7_KO_Zuzka/191125_A00380_0034_AHHT2MDRXX"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- file.path(inpath, "190925_ZL_PAPD7 testes and liver libraries preparation.xlsx")

######################################################## READ DATA
# read raw sample data
sample_table_raw <- openxlsx::read.xlsx(sample_table_raw_path)

######################################################## MAIN CODE
# get experiment name
experiment <- "Papd7_KO"
  
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name = X3, tissue = X6, genotype = X4, barcode = X15, bc = X14) %>% 
  dplyr::slice(-1) %>% 
  dplyr::filter(!is.na(sample_name)) %>% 
  dplyr::mutate(sample_name = str_replace_all(sample_name, "-", "_") %>% str_trim(.) %>% str_c(., str_extract(tissue, ".{1}"))) 

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T, recursive = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz"), 
                sample_name = str_extract(sample_file_name, pattern = str_c(sample_table_tidy$sample_name, collapse = "|")))

# check if you have all samples
if(!(all(sample_tb$sample_name %in% sample_table_tidy$sample_name))){
  
  # warn me
  warning("You lost some samples mate! Check what's wrong!")
  
}

# join by sample name
sample_tb %<>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_name") %>% 
  dplyr::mutate(sample_name = str_remove(sample_name, "[l,t]$"), 
                genotype = str_replace_all(genotype, " ", "_") %>% str_remove(., "\\?")) %>% 
  dplyr::group_by(tissue, genotype) %>% 
  dplyr::mutate(replicate = 1:n(), 
                sample_id = str_c("s", tissue, genotype, sample_name, str_c("r", replicate), sep = "_") %>% str_c(., ".SE"), 
                sample_id = replace(sample_id, sample_file_name == "Undetermined_S0_R1_001", "s_undetermined_r1.SE"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, tissue, genotype, replicate, bc, barcode, sample_name, sample_file_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = "."))) %T>% 
  readr::write_csv(., path = file.path(raw_path, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))


