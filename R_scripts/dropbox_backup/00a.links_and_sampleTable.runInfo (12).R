### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Wang_2013_EmboJ_GSE41455"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/ChIPseq", reference, "Data/Raw/Links"))

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
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/ChIPseq", reference)

# documentation path
runinfo_path <- list.files(path = sra_path, pattern = "*runInfo.txt", full.names = T)

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = runinfo_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# geo table
geo_tb <- tibble(geo = c("GSM1017630", 
                         "GSM1017632", 
                         "GSM1017634", 
                         "GSM1017635"), 
                 sample_type = c("H3K4me3_ChIPSeq", 
                                 "H3K9me2_ChIPSeq", 
                                 "Input for H3K4me3 and H3K27me3", 
                                 "Input for H3K9me2"))

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run, 
    geo = `GEO_Accession (exp)`, 
    antibody = Antibody
    
  ) %>% 
  
  dplyr::inner_join(., geo_tb, by = "geo") %>% 
  # dplyr::filter(!str_detect(sample_type, "\\(d7\\)")) %>% 
  dplyr::mutate(sample_type = str_remove_all(sample_type, " for| and"), 
                sample_type = str_replace(sample_type, "Input", "input"), 
                sample_type = str_replace_all(sample_type, " ", "_"), 
                antibody = ifelse(antibody == "none", sample_type, antibody), 
                antibody = replace(antibody, antibody == "input_H3K4me3_H3K27me3", "input_H3K4me3")) %>% 
  dplyr::select(-sample_type) %>% 
  
  dplyr::mutate(
    
    cell = "3T3"
    
  ) %>% 
  
  dplyr::mutate(name = str_c("s", 
                             cell,
                             antibody,
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
                cell,
                antibody) %T>%
  readr:::write_csv(., file = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., file = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

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
    arrange(run)
  
}


# save make links script
rename_file %<>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., file = file.path(outpath, "make_links.sh"))
