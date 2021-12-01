### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Gahurova_2017_EpigeneticsChromatin_GSE86297"

# set working directory
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/otherSeq", reference, "Data/Raw/Links"))

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
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/otherSeq", reference)

# documentation path
runinfo_path <- list.files(path = sra_path, pattern = "*runInfo.txt", full.names = T)

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = runinfo_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
# pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run, 
    stage = Developmental_stage, 
    size = source_name, 
    pairing = LibraryLayout
    # genotype = Genotype, 
    # antibody = Antibody, 
    # assay = `Assay Type`
    # age, 
    # genotype = strain
    
  ) %>% 
  
  dplyr::mutate(
    
    stage = "oocyte", 
    size = str_remove(size, " oocytes"), 
    pairing = str_replace_all(pairing, c("SINGLE" = "SE", "PAIRED" = "PE")), 
    library_type = ifelse(pairing == "SE", "RRBS", "PBAT") 

  ) %>% 
  

  dplyr::mutate(name = str_c("s", 
                             stage,
                             size,
                             library_type,
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
                size, 
                pairing, 
                library_type, 
                run) %T>%
  readr:::write_csv(., path = file.path(outpath, "../../Documentation", str_c(reference, ".sampleTable.csv"))) %T>%
  readr:::write_csv(., path = file.path(sra_path, str_c(reference, ".sampleTable.csv")))

# only for PE data
if(pairing == "PE"){

  rename_file %>%
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
rbind(rename_file,
      rename_file %>% dplyr::filter(pairing == "PE")) %>% 
  dplyr::mutate(name = ifelse(pairing == "PE", str_c(name, "_", 1:2), name)) %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = str_c(outpath, "/", name, ".txt.gz"),
                run = ifelse(pairing == "PE", str_c(sra_path, "/", run, "_", 1:2, ".fastq.gz"), str_c(sra_path, "/", run, ".fastq.gz"))) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", run, " ", name)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
