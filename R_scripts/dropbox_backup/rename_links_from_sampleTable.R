### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/golden_hamster.oocytes_testis.2019/190923_A00380_0028_BHFY7GDRXX_S016_malik"

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
  basename(.) %>% 
  str_remove(., "\\.RNAseq|\\.small_RNAseq")
  
# clean sample table
sample_table_tidy <- 
  sample_table_raw %>% 
  dplyr::select(sample_name = sample, genotype = `genotype Mov10l`, barcode = `Barcode seq`, bc = `Barcode number`)

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz"), 
                sample_name = str_extract(sample_file_name, pattern = str_c(sample_table_tidy$sample_name, collapse = "|"))) %>% 
  dplyr::filter(!is.na(sample_name)) %>% 
  dplyr::filter(!str_detect(sample_file_name, "-M[0-9]+"))

# check if you have all samples
if(!(all(sample_tb$sample_name %in% sample_table_tidy$sample_name))){
  
  # warn me
  warning("You lost some samples mate! Check what's wrong!")
  
}

# join by sample name
sample_tb %<>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_name") %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_c("Mov10l_", genotype)) %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::mutate(replicate = 1:n(), 
                sample_id = str_c("s", stage, genotype, str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, bc, barcode, sample_name, sample_file_name) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



