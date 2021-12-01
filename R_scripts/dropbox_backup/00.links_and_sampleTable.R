### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Raw/Links")

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

# raw path
raw_path <- "/common/RAW/Svoboda/golden_hamster.bisulfite_seq.test_run/2020_Sep"

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Documentation"
sample_table_path <- list.files(path = documentation_path, pattern = "*sampleTable.raw.csv", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "PE"

# create rename executable
rename_file <- 
  
  sample_table %>% 
  
  dplyr::select(
    
    sample_name = sample, 
    species, 
    number_of_oocytes = oocytes,
    number_of_PCR_cycles = `number of PCR cycles`, 
    size_selection = `size selection`, 
    barcode = `Barcode seq`
    
  ) %>% 
  
  dplyr::mutate(
    
    stage = "GV", 
    genotype = ifelse(species == "mouse", "WT", "Mov10l1_HET"), 
    species = str_remove(species, " .*")

  ) %>% 
  

  dplyr::mutate(name = str_c("s", 
                             species, 
                             stage,
                             genotype, 
                             sample_name, 
                             sep = "_")) %>% 
  
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %<>%
  dplyr::select(sample_id = name, 
                species, 
                stage,
                genotype, 
                sample_name, 
                barcode, 
                everything()) %T>%
  readr:::write_csv(., path = sample_table_path %>% str_replace_all(., " |-", "_") %>% str_replace(., "\\.sampleTable\\.raw\\.csv$", ".sampleTable.csv"))

# join with raw file names
raw_file_names <- 
  tibble(raw_path = list.files(raw_path, "\\.fastq\\.gz", full.names = T, recursive = F)) %>% 
  dplyr::mutate(species = str_extract(raw_path, "Hamster|Mouse") %>% tolower(.), 
                read_in_pair = str_extract(raw_path, "(?<=r)[1,2]{1}(?=\\.fastq\\.gz)") %>% as.integer(.))

# only for PE data
if(pairing == "PE"){

  rename_file %<>%
    rbind(., .) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(read_in_pair = 1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., raw_file_names, by = c("species", "read_in_pair")) %>% 
    dplyr::mutate(name = str_c(sample_id, "_", read_in_pair)) %>%
    dplyr::select(name, raw_path) %>%
    arrange(name)

}


# save make links script
rename_file %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
