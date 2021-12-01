### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Jin_2021_NatCommun_GSE162148"

# set working directory
setwd(file.path("."))

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
runinfo <- readr::read_csv(file = runinfo_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run,
    tissue, 
    age = AGE, 
    sex,
    organism = Organism, 
    breed, 
    assay = `Assay Type`, 
    library_name = `Library Name`
    
  ) %>% 
  
  dplyr::filter(organism == "Homo sapiens") %>%
  dplyr::group_by(organism) %>% 
  dplyr::count(tissue)
  
  dplyr::mutate(tissue = str_replace_all(tissue, " ", "_") %>% tolower(.), 
                experiment = "sRNAseq") %>% 
  
  dplyr::mutate(name = str_c("s", 
                             tissue,
                             experiment,
                             sep = "_")) %>% 
  dplyr::group_by(name) %>% 
  dplyr::mutate(sequence = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>%
  dplyr::select(sample_id = name, 
                tissue,
                experiment, 
                run) %T>%
  readr:::write_csv(., file = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., file = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

# # only for PE data
# if(pairing == "PE"){
#   
#   rename_file %<>% 
#     rbind(., .) %>% 
#     dplyr::group_by(run) %>% 
#     dplyr::mutate(sequence = 1:n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(run = str_c(run, "_", sequence), 
#                   name = str_c(name, "_", sequence)) %>% 
#     dplyr::select(-sequence) %>% 
#     arrange(run)
#   
# }


# save make links script
rename_file %<>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., file = file.path(outpath, "make_links.sh"))
