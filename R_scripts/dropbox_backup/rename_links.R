### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "CNOT6L"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Raw/Links"))

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
sra_path <- "/common/RAW/Svoboda/HamsterUtah/11919R/Fastq"

# documentation path
runinfo_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Documentation"

# sample table path
sample_table_path <- file.path(runinfo_path, "CNOT6L_sample_list_11919R_2015_10_29.csv")

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = sample_table_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "PE"

# raw files
raw_table <- 
  list.files(sra_path, ".*\\.txt\\.gz$", full.names = T) %>% 
  tibble(raw_path = .) %>% 
  dplyr::mutate(file_name = raw_path %>% basename(.),
                run = file_name %>% str_remove(., "_.*$")) %>% 
  dplyr::mutate(run = str_c(run, "_", 1:2)) %>%  
  dplyr::select(run, file_name, raw_path)

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = ID, 
    stage = `Time Course`,
    genotype = `Treatment/Control`
    # genotype = Genotype, 
    # antibody = Antibody, 
    # assay = `Assay Type`
    # age, 
    # genotype = strain
    
  ) %>% 
  
  dplyr::mutate(stage = factor(stage, levels = c("GV", "MII", "1C")), 
                genotype = factor(genotype, levels =c("WT", "KO", "Hamster"))) %>% 

  dplyr::arrange(stage, genotype) %>% 
  
  dplyr::mutate(name = str_c("s", 
                             stage,
                             genotype,
                             sep = "_")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>%
  dplyr::select(sample_id = name, 
                stage,
                # cell,
                genotype,
                # antibody, 
                run) %T>%
  readr:::write_csv(., path = file.path(runinfo_path, str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., path = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

# only for PE data
if(pairing == "PE"){
  
  rename_file %<>% 
    rbind(., .) %>% 
    dplyr::group_by(run) %>% 
    dplyr::mutate(sequence = 1:n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(run = str_c(run, "_", sequence), 
                  name = str_c(name, "_", sequence)) %>% 
    dplyr::select(-sequence) %>% 
    dplyr::left_join(., raw_table, by = "run")

}


# save make links script
rename_file %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
