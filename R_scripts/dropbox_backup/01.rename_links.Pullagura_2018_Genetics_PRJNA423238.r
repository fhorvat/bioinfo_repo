### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Pullagura_2018_Genetics_PRJNA423238"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq", reference, "Data/Raw/Links"))

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
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/smallRNAseq", reference)

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
    sample_name = `Library Name`, 
    sample_name_short = `Sample Name`
    # genotype = Genotype, 
    # antibody = Antibody, 
    # assay = `Assay Type`
    # age, 
    # genotype = strain
    
  ) %>% 
  
  dplyr::mutate(
    
    genotype = str_extract(sample_name, "Tarbp2_Mut|Prkra_Mut|Tarbp2_WT|Prkra_WT"), 
    replicate = str_extract(sample_name, "[0-9]+$") %>% as.integer(.)
    
  ) %>%
  
  dplyr::arrange(genotype, replicate) %>% 

  dplyr::mutate(
    
    stage = "embryo"
    
  ) %>% 
  
  # dplyr::arrange(run) %>% 
  
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
                sample_name, 
                sample_name_short,
                run) %T>%
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
