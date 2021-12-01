### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_mESC.LCMV.small_RNAseq.2021_Sep/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/mouse_mESC.LCMV.small_RNAseq.2021_Sep"

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
  dplyr::select(sample_name = sample, genotype, 
                experiment = Experiment, 
                animal_or_cell_ID = `ID (cell samples)`,
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  dplyr::mutate(experiment_short = str_remove(experiment, ", .*") %>% str_replace(., " ", "_"),  
                infection = "LCMV",
                genotype = str_replace(genotype, "Dcr", "Dicer"),
                genotype = str_replace(genotype, "DicerOO, PKR del", "DicerOO_PKRdel"), 
                animal_or_cell_ID = as.character(animal_or_cell_ID)) %>% 
  dplyr::group_by(sample_name, experiment_short) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", sample_name, animal_or_cell_ID, experiment_short, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_name = basename(raw_path)) %>% 
  dplyr::filter(str_detect(sample_name, "^RS")) %>%
  dplyr::mutate(sample_name = str_remove(sample_name, "\\.fastq\\.gz") %>% 
                  str_remove(., "_S[0-9]+$"), 
                animal_or_cell_ID = str_extract(sample_name, "(?<=LCMV)[0-9]+$"), 
                sample_name = str_remove(sample_name, "(?<=LCMV)[0-9]+$") %>% 
                  str_replace_all(., "-", "_")) %>% 
  dplyr::left_join(., sample_table_tidy, by = c("sample_name", "animal_or_cell_ID"))


# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, genotype, infection, experiment, replicate, 
                sample_name, animal_or_cell_ID,
                barcode, barcode_no, raw_path) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



