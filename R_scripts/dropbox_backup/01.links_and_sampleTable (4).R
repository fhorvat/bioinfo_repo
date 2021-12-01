### INFO: 
### DATE: Wed Nov 27 15:10:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Nov/Data/Documentation")

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

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# raw files path
raw_path <- "/common/RAW/Svoboda/DicerX_viral_infection/demultiplexed_reads"

# links path 
links_path <- file.path(inpath, "../Raw/Links") 

# sample table path 1
sample_table_raw_path_1 <- file.path(inpath, "DicerX_spleen_smallRNA_libraries_seqID_191122.sheet1.csv")

# sample table path 2
sample_table_raw_path_2 <- file.path(inpath, "DicerX_spleen_smallRNA_libraries_seqID_191122.sheet2.csv")


######################################################## READ DATA
# read raw sample data 1
sample_table_raw_1 <- readr::read_csv(sample_table_raw_path_1)

# read raw sample data 2
sample_table_raw_2 <- readr::read_csv(sample_table_raw_path_2)

######################################################## MAIN CODE
# clean sample table 1
sample_table_tidy_1 <- 
  sample_table_raw_1 %>% 
  dplyr::rename_all(str_replace_all, " ", "_") %>% 
  dplyr::select(sample_name = Sample_ID, genotype = Genotype, treatment = Treatment, barcode) %>% 
  dplyr::mutate(genotype = str_extract(genotype, "WT|HET"), 
                barcode = as.character(barcode))

# clean sample table 2
sample_table_tidy_2 <- 
  sample_table_raw_2 %>% 
  dplyr::rename_all(str_replace_all, " ", "_") %>% 
  dplyr::select(sample_file_name = Sample_ID, barcode = Barcode_ID, barcode_seq = Barcode_seq) %>% 
  dplyr::mutate(barcode = str_remove(barcode, "PCR Primer "), 
                barcode_seq_rc = DNAStringSet(barcode_seq) %>% reverseComplement(.) %>% as.character(.))

# create sample table
sample_tb <- 
  tibble(raw_path = list.files(raw_path, pattern = "*.\\.txt\\.gz", full.names = T)) %>% 
  dplyr::mutate(barcode_seq_rc = raw_path %>% basename(.) %>% str_remove_all(., "out_|\\.txt\\.gz")) %>% 
  dplyr::left_join(., sample_table_tidy_2, by = "barcode_seq_rc") %>% 
  dplyr::left_join(., sample_table_tidy_1, by = "barcode") %>% 
  dplyr::mutate(stage = "spleen") %>% 
  dplyr::group_by(genotype, treatment, stage) %>% 
  dplyr::mutate(replicate = 1:n()) %>%
  dplyr::ungroup(.) %>% 
  dplyr::select(sample_name, sample_file_name, genotype, treatment, stage, replicate, barcode_seq, barcode_seq_rc, barcode, raw_path) %>% 
  dplyr::arrange(genotype, treatment, stage, replicate) %>% 
  dplyr::mutate(sample_id = str_c("s", stage, "DicerX", genotype, treatment, sample_name, sep = "_") %>% str_c(., "_r", replicate, ".SE"), 
                sample_id = replace(sample_id, is.na(sample_id), "s_spleen_unmatched_r1.SE"))

# save as renaming script
sample_tb %>%
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", links_path, "/", sample_id, ".txt.gz")) %$% 
  make_links %T>%
  readr::write_lines(., path = file.path(links_path, "make_links.sh"))

# save as sample table
sample_tb %>% 
  dplyr::select(sample_id, stage, genotype, treatment, replicate, sample_name, barcode, barcode_seq, barcode_seq_rc, sample_file_name, raw_path) %>% 
  dplyr::distinct(.) %>% 
  dplyr::arrange(sample_id) %T>% 
  readr:::write_csv(., path = file.path(outpath, str_c("DicerX_spleen_smallRNA_libraries_seqID_191122", ".sampleTable.csv")))







