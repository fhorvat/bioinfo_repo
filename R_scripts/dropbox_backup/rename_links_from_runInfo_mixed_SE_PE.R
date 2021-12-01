### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Kabayama_2017_NucleicAcidsRes_PRJDB4628"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", reference, "Data/Raw/Links"))

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
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/RNAseq", reference)

# documentation path
runinfo_path <- list.files(path = sra_path, pattern = "*runInfo.txt", full.names = T)

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_delim(file = runinfo_path, delim = "\t")

######################################################## MAIN CODE
# # set whether reads are single or paired end
# pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    run = Run, 
    sample_name = Library_Name,
    tissue = tissue_type, 
    layout = LibraryLayout
  ) %>% 
  
  dplyr::mutate(
    genotype = 
      str_extract_all(sample_name, "Mili KO|WT") %>% 
      str_replace_all(., " ", "_"), 
    layout = str_replace_all(layout, c("SINGLE" = "SE", "PAIRED" = "PE"))
  ) %>%
  
  dplyr::mutate(name = str_c("s", 
                             # cell, 
                             # tissue,
                             genotype,
                             sep = "_")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", layout)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>%
  dplyr::select(sample_id = name, 
                tissue,
                # stage,
                genotype,
                run) %T>%
  readr:::write_csv(., path = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., path = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

# PE and SE data
rename_file <- 
  
  rbind(
  
  rename_file %>% 
    dplyr::filter(layout == "PE") %>% 
    rbind(., .) %>% 
    dplyr::group_by(run) %>% 
    dplyr::mutate(sequence = 1:n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(run = str_c(run, "_", sequence), 
                  name = str_c(name, "_", sequence)) %>% 
    dplyr::select(-sequence) %>% 
    arrange(run), 
  
  rename_file %>% 
    dplyr::filter(layout == "SE")
  
)

# save make links script
rename_file %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))

