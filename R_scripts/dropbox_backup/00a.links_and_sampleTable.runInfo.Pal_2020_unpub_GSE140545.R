### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Pal_2020_unpub_GSE140545"

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

# GSM4173984 	Sample 1: Induced pluripotent stem cell
# GSM4173985 	Sample 2: Retinal progenitors
# GSM4173986 	Sample 3: RPE progenitors
# GSM4173987 	Sample 4: Retinal Pigment Epithelium
# GSM4173988 	Sample 5: Photoreceptor precursors

# create geo table
geo_tb <- tibble(geo_sample = c("GSM4173984", "GSM4173985", "GSM4173986", "GSM4173987", "GSM4173988"), 
                 sample = c("Induced pluripotent stem cell", "Retinal progenitors", "RPE progenitors", 
                            "Retinal Pigment Epithelium", "Photoreceptor precursors"), 
                 sample_name = c("iPSC", "retinal_progenitors", "RPE_progenitors",
                                 "RPE", "photoreceptor_precursors"))

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run,
    geo_sample = `Sample Name`
    
  ) %>% 
  
  dplyr::left_join(., geo_tb, by = "geo_sample") %>% 

  dplyr::mutate(name = str_c("s", 
                             sample_name,
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
                sample_name,
                sample, 
                geo_sample,
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
