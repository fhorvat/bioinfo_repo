### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/rodent_oocytes.small_RNAseq.2021_Sep/Data/Documentation")

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

# raw files path
raw_path <- "/common/RAW/Svoboda/rodent_oocytes.small_RNAseq.2021_Sep"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw\\.csv", full.names = T)

######################################################## READ DATA
# read raw sample data
sample_table_raw <- readr::read_csv(sample_table_raw_path)

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.)

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::mutate(genotype = str_extract(strain, "WT|Dicer.*|Mov10l1.*"), 
                genotype = str_remove(genotype, " hamster"),
                genotype = replace(genotype, is.na(genotype), "WT"), 
                genotype = str_replace_all(genotype, " ", "_") %>% 
                  str_replace(., "DicerSOM", "Dicer_SOM") %>% 
                  str_replace(., "\\+/\\+", "Homozyg") %>% 
                  str_replace(., "-/-", "KO")) %>% 
  dplyr::select(sample, tissue, genotype, strain, age = `age (dpp)`, barcode = `Barcode seq`) %>% 
  dplyr::mutate(age = tolower(age),
                strain = str_remove(strain, " WT"), 
                strain = replace(strain, strain == "BL6", "BL6Crl"),
                strain = replace(strain, str_detect(strain, "Dicer"), "BL6Crl"), 
                strain = replace(strain, str_detect(strain, "Mov10l1"), "hamster"), 
                strain = replace(strain, strain == "Castaneus", "Cast"),
                strain = str_remove(strain, ":"), 
                sample_name = str_remove_all(sample, ".*_| ham KO .*")) %>% 
  dplyr::group_by(strain, genotype) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = ifelse(strain != "hamster", 
                                   str_c("s", strain, genotype, sep = "_"),
                                   str_c("s", genotype, sample_name, sep = "_")), 
                sample_id = str_c(sample_id, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = ".*\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz")) %>% 
  dplyr::mutate(sample_name = ifelse(str_detect(sample_file_name, "^So[0-9]{3}"),
                                     str_remove(sample_file_name, "-.*oocytes"), 
                                     str_remove(sample_file_name, ".*-"))) %>% 
  dplyr::mutate(sample_name = str_replace(sample_name, "04734", "812-04734")) %>% 
  dplyr::mutate(sample_name = str_remove(sample_name, "_.*S[0-9]+")) %>% 
  dplyr::left_join(., sample_table_tidy, by = c("sample_name"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, strain, genotype, age, tissue, replicate, sample, sample_name, barcode, raw_path) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



