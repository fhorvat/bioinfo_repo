### INFO: 
### DATE: Sun Sep 22 23:11:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc5_KO.2020_Dec/Data/Documentation")

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
raw_path <- "/common/RAW/Svoboda/lncRNA_KO/201217_Linc5del"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path
sample_table_path <- list.files(inpath, ".*\\.sampleTable\\.raw\\.csv")

######################################################## READ DATA
# read sample table
sample_table_raw <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
# get experiment name
experiment <- 
  outpath %>% 
  str_remove(., "/Data/.*$") %>% 
  basename(.) 

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.fastq\\.gz", full.names = T)) %>% 
  dplyr::mutate(file_name = raw_path %>% basename(.) %>% str_remove_all(., "\\.fastq\\.gz"), 
                file_name = str_replace_all(file_name, "-", "_") %>% str_remove(., "_S[0-9]{1}")) %>% 
  dplyr::left_join(., sample_table_raw, by = c("file_name" = "Name")) %>% 
  dplyr::mutate(stage = "GV",
                genotype = str_extract(file_name, "Bl6_wt|Linc5_del|Linc5_wt") %>% toupper(.)  %>% 
                  str_replace(., "LINC", "Linc"), 
                sample_replicate = SampleID %>% str_remove(., ".*_") %>% str_c("S", .)) %>% 
  dplyr::group_by(stage, genotype) %>%
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sample_id = str_c("s", stage, genotype, sample_replicate, str_c("r", replicate), sep = "_") %>% str_c(., ".SE"))

# save as renaming script
links <- 
  sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., file = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, replicate, file_name, Index1Name, Index1Sequence, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., file = file.path(outpath, str_c(experiment, 
                                                       format(Sys.Date(), "%Y%m%d"), 
                                                       "sampleTable.csv", sep = ".")))



