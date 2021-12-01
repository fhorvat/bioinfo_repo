### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "Li_2017_MolCell_GSE90874"

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
runinfo_path <- list.files(path = sra_path, pattern = "*\\.runInfo\\.txt", full.names = T)
geoinfo_path <- list.files(path = sra_path, pattern = "*\\.geoInfo\\.txt", full.names = T)

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = runinfo_path)
geoinfo <- readr::read_delim(file = geoinfo_path, delim = "\t", col_names = c("geo", "sample_name"))

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- ifelse((runinfo$LibraryLayout %>% unique(.)) == "SINGLE", "SE", "PE")

# create rename executable
rename_file <- 
  
  runinfo %>% 
  
  dplyr::select(
    
    run = Run, 
    source = source_name, 
    treatment,
    type, 
    geo = `GEO_Accession (exp)`
    # genotype = Genotype, 
    # antibody = Antibody, 
    # assay = `Assay Type`
    # age, 
    # genotype = strain
    
  ) %>% 
  
  dplyr::mutate(
    
    source = "HeLa", 
    type = replace(type, (is.na(type) | type == "HeLa with expressing Flag-NF90"), "Flag_NF90_iCLIP"), 
    type = str_replace_all(type, c("poly\\(A\\)\\+ RNA" = "polyA_pos_RNA", 
                                   "steady-state ribo-\\/poly\\(A\\)- RNA" = "ribo_neg_polyA_neg_RNA", 
                                   "nascent RNA labeled with 4sU 2hr ribo- RNA" = "nascent_4sU_ribo_neg_RNA"))
    
  ) %>%
  
  dplyr::left_join(., geoinfo, by = "geo") %>% 
  
  dplyr::mutate(genotype = str_extract(sample_name, "NF90 KD|NF110 KD|scramble|scramble control") %>% 
                  str_replace_all(., " ", "_") %>% 
                  str_replace_all(., "control", "ctrl")) %>% 
  
  dplyr::mutate(genotype = replace(genotype, is.na(genotype), "")) %>% 
  
  dplyr::arrange(run) %>% 
  
  dplyr::mutate(name = str_c("s", 
                             source,
                             genotype, 
                             type,
                             sep = "_")) %>% 
  dplyr::mutate(name = str_replace(name, "__", "_")) %>% 
  dplyr::group_by(name) %>%
  dplyr::mutate(sequence = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = str_c(name, "_r", sequence),
                name = str_c(name, ".", pairing)) %>%
  dplyr::select(-sequence)

# save sample table
rename_file %>%
  dplyr::select(sample_id = name, 
                sample_name,
                genotype, 
                type,
                treatment, 
                source,
                # cell,
                # genotype,
                # antibody, 
                geo,
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
