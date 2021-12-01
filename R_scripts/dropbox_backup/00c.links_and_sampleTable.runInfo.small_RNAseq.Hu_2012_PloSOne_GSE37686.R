### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Hu_2012_PloSOne_GSE37686"

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
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/small_RNAseq", reference)

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
    cell_type = Cell_type
    
  ) %>% 
  
  dplyr::mutate(sample_name = str_remove(sample_name, ".*: ") %>% 
                  str_replace(., "partial differentiated|partially differentiated", "part_diff") %>% 
                  str_replace_all(., " ", "_") %>% 
                  str_replace(., "ES-RPE", "ES_RPE") %>% 
                  str_replace(., "part_diff_ESC_BG01_partial_diff", "ESC_BG01_part_diff") %>% 
                  str_remove(., "_rep[1,2]") %>% 
                  str_replace(., "part_diff_ESC_H9", "ESC_H9_part_diff") %>% 
                  str_replace(., "fetal_RPE", "prim_fetal_RPE") %>% 
                  str_replace(., "pigmented_cluster", "pigment_clust")) %>% 
  
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
                cell_type, 
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
