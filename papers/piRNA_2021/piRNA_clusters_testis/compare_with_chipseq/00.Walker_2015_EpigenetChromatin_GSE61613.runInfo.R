### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Walker_2015_EpigenetChromatin_GSE61613"

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
geoinfo_path <- list.files(path = sra_path, pattern = "*geoInfo.txt", full.names = T)

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = runinfo_path)
geoinfo <- readr::read_csv(file = geoinfo_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# tidy geo info
geoinfo_tb <- 
  geoinfo %>%
  tidyr::separate(`Experiment Title`, into = c("geo", "sample", "tmp1", "tmp2"), sep = ";|:") %>% 
  dplyr::mutate_all(~str_trim(.)) %>% 
  dplyr::select(1:2)

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::filter(`Assay Type` == "ChIP-Seq") %>% 
  
  dplyr::select(
    
    run = Run, 
    geo = `GEO_Accession (exp)`, 
    age = Age, 
    antibody = Antibody, 
    tissue = source_name
    
    
  ) %>% 
  
  dplyr::left_join(., geoinfo_tb, by = "geo") %>% 
  dplyr::mutate(age = str_replace(age, " days post-partum", "dpp"), 
                sample = replace(sample, str_detect(sample, "GES14"), "H3K9me2_ChIPSeq") %>% 
                  str_replace(., "Input DNA", "input") %>% 
                  str_replace(., "Affinity-seq", "affinity") %>% 
                  str_replace(., "ChIPSeq", "ChIP") %>% 
                  str_replace_all(., " ", "_") %>% 
                  replace(., . == "input", "H3K9me2_input") %>% 
                  replace(., . == "affinity_input", "Prdm9_affinity_input") %>% 
                  replace(., . == "Prdm9_input", "Prdm9_ChIP_input")) %>% 
  dplyr::mutate(name = str_c("s", 
                             sample,
                             age,
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
                name = sample, 
                antibody,
                tissue,
                age,
                run, 
                geo) %T>%
  readr:::write_csv(., path = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
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
    arrange(run)
  
}


# save make links script
rename_file %<>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = str_c(sra_path, "/", run, ".fastq.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
