### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.mESC.TBEV_LCMV.small_RNAseq.2021_May/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse_brain.mESC.TBEV_LCMV.small_RNAseq.2021_May"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw\\.csv", full.names = T)

######################################################## READ DATA
# read raw sample data
sample_table_raw <- data.table::fread(sample_table_raw_path)

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.)

# clean sample table
sample_table_tidy <- 
  sample_table_raw %>%
  as_tibble(.) %>% 
  dplyr::filter(sample != "") %>%
  dplyr::select(sample_name = sample, genotype, animal_or_cell_ID = `ID (animals or cell samples)`,
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  dplyr::mutate(tissue = str_extract(sample_name, "^RS[0-9]+"), 
                tissue = ifelse(is.na(tissue), "brain", str_c(tissue, "_mESC")),
                infection = str_extract(sample_name, "TBEV|LCMV"),
                genotype = str_replace(genotype, "Dcr", "Dicer"),
                genotype = str_replace(genotype, "DicerOO, PKR del", "DicerOO_PKRdel"), 
                animal_or_cell_ID = as.character(animal_or_cell_ID)) %>% 
  dplyr::group_by(tissue, genotype, infection) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", tissue, genotype, infection, animal_or_cell_ID, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = basename(raw_path), 
                genotype = str_extract(sample_file_name, "DicerX|DcrOO-PKRdel|WT") %>% 
                  str_replace(., "Dcr", "Dicer") %>% 
                  str_replace(., "-", "_"),
                infection = str_extract(sample_file_name, "TBEV[0-9]+|LCMV[0-9]+"), 
                animal_or_cell_ID = str_extract(infection, "[0-9]+"), 
                infection = str_remove(infection, "[0-9]+"), 
                tissue = str_extract(sample_file_name, "^RS[0-9]+"), 
                tissue = ifelse(is.na(tissue), "brain", str_c(tissue, "_mESC"))) %>% 
  dplyr::left_join(., sample_table_tidy, by = c("tissue", "genotype", "infection", "animal_or_cell_ID"))


# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, tissue, genotype, infection, replicate, 
                sample_name, sample_file_name, animal_or_cell_ID,
                barcode, barcode_no, raw_path) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



