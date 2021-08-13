### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Hendrickson_2017_NatGenet_GSE72379"

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
runinfo_path <- 
  list.files(path = sra_path, pattern = "*.runInfo.txt", full.names = T) %>% 
  .[!str_detect(., "GV")]

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_delim(file = runinfo_path, delim = "\t")

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run, 
    stage = developmental_stage, 
    name = combined_analysis_name

  ) %>% 
  
  dplyr::mutate(stage = str_replace_all(stage, c("1-cell.*" = "1C", 
                                                 "2-,4-,.*" = "2C_4C_8C", 
                                                 "Blastocyst; ICM.*" = "blast_ICM_polarTroph", 
                                                 "Blastocyst; mural.*" = "blast_muralTroph", 
                                                 "Morula" = "morula")), 
                stage = ifelse(stage %in% c("Immature oocyte", "Mature oocyte"), name, stage), 
                stage = str_remove(stage, "_rep.*"), 
                technical_replicate = str_extract(name, "rep[0-9]{1}$") %>% str_replace(., "rep", "tr")) %>% 
  
  dplyr::mutate(name = str_c("s", 
                             stage,
                             technical_replicate,
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
                technical_replicate,
                run) %T>%
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
