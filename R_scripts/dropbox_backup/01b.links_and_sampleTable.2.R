### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.0.5dpp.RNAseq/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/hamster_testis_ovaries_Mov10l.0.5dpp.smallRNAseq"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_raw_path <- list.files(inpath, pattern = ".*\\.sampleTable\\.raw_2\\.csv", full.names = T)

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
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::select(sample, tissue, genotype = `genotype Mov10l`, age = `age (dpp)`,
                barcode = `Barcode seq`, barcode_no = `Barcode number`) %>% 
  dplyr::mutate(tissue = str_replace(tissue, "testes", "testis") %>% str_remove(., " "),
                tissue = replace(tissue, str_detect(tissue, "totalRNA"), "testis_RNAseq"),  
                tissue = str_replace(tissue, "oocyes", "oocytes"), 
                age = str_replace(age, ",", "."),
                sample_name = str_remove_all(sample, "\\(|\\)") %>% 
                  str_replace(., " ", "_") %>% 
                  str_replace(., "SO(?=[0-9]+)", "So"), 
                genotype = str_c("Mov10l1_", genotype)) %>% 
  dplyr::group_by(genotype, tissue) %>%
  dplyr::mutate(replicate = 1:n(),
                sample_id = str_c("s", tissue, genotype, age, sample_name, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = ".*\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(sample_file_name = raw_path %>% basename(.) %>% str_remove_all(., "_r1\\.fastq\\.gz"), 
                # genotype = str_extract(sample_file_name, "WT|HET|KO") %>% str_c("Mov10l1_", .), 
                # genotype = replace(genotype, is.na(genotype), "Mov10l1_WT"),
                tissue = str_extract(sample_file_name, "[0-9]+oocytes|ovary"),
                tissue = replace(tissue, is.na(tissue), "testis"), 
                tissue = replace(tissue, str_detect(sample_file_name, "total"), "testis_RNAseq")) %>% 
  dplyr::mutate(sample_name = str_remove(sample_file_name, "-WT.*|-HET.*|-KO.*") %>% 
                  str_remove(., "_.*") %>% 
                  str_replace(., "oocytes", "xoo") %>% 
                  str_replace(., "-(?=[0-9]+)", "_") %>% 
                  str_remove(., "-ovary") %>% 
                  str_replace(., "19xoo", "19x")) %>% 
  dplyr::left_join(., sample_table_tidy, by = c("sample_name", "tissue"))

# filter
sample_tb %<>% 
  dplyr::filter(tissue == "testis_RNAseq") %>% 
  dplyr::mutate(tissue = "testis", 
                replicate = 2, 
                genotype = str_replace(genotype, "Mov10l1_", "Mov10l_")) %>% 
  dplyr::mutate(sample_id = str_c("s", tissue, genotype, age, sample_name, 
                                  str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))


# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.2.sh"))

# save as sample table
sample_tb %<>% 
  dplyr::select(sample_id, tissue, genotype, age, replicate, barcode, sample_name = sample) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.part_2.csv", sep = ".")))

# read part 1 of sample table, join all together, save
readr::read_csv(file.path(inpath, "hamster_testis_Mov10l.0.5dpp.RNAseq.20210428.sampleTable.part_1.csv")) %>% 
  dplyr::bind_rows(., sample_tb) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))


