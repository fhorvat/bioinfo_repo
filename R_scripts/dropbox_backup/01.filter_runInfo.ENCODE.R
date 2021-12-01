### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
# set data/reference
reference <- "ENCODE_2014_Nature_GSE49417"

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
# inpath
inpath <- getwd()

# outpath
outpath <- getwd()

# SRA path
sra_path <- file.path("/common/DB/SRA/sra/Svoboda/2017_download/RNAseq", reference)

# documentation path
runinfo_path <- list.files(path = sra_path, pattern = "ENCODE_2014_Nature_GSE49417\\.runInfo\\.txt", full.names = T)

# sample table path
sample_table_path <- file.path(inpath, "../../Documentation")
sample_table_path <- list.files(sample_table_path, "\\.sampleTable.csv", full.names = T)

# links path
links_path <- file.path(inpath, "make_links.sh")

######################################################## READ DATA
# read runInfo.txt table
runinfo <- readr::read_csv(file = runinfo_path)

# read sample table
sample_tb <- readr::read_csv(sample_table_path)
  
# links lines
links_tb <- readr::read_delim(links_path, delim = " ", col_names = c("tmp_1", "tmp_2", "sra", "sample_id"))

######################################################## MAIN CODE
# create sample id list
sample_id_list <- c("s_liver_adult8wks_r1.PE",
                    "s_liver_adult8wks_r2.PE",
                    "s_liver_adult8wks_r3.PE",
                    "s_liver_adult8wks_r4.PE",
                    "s_liver_adult8wks_r5.PE", 
                    "s_liver_adult8wks_r6.PE", 
                    "s_heart_adult8wks_r1.PE", 
                    "s_heart_adult8wks_r2.PE", 
                    "s_heart_adult8wks_r3.PE", 
                    "s_heart_adult8wks_r4.PE")

# clean links table
links_tb_tidy <- 
  links_tb %>%
  dplyr::select(-c(tmp_1, tmp_2)) %>% 
  dplyr::mutate(sra = basename(sra) %>% str_remove(., "_[1,2]{1}\\.fastq\\.gz"),
                sample_id = basename(sample_id) %>% str_remove(., "_[1,2]{1}\\.txt\\.gz")) %>% 
  unique(.) %>% 
  dplyr::filter(sample_id %in% sample_id_list)

# filter runInfo
runinfo_filt <- 
  runinfo %>% 
  dplyr::filter(Run %in% links_tb_tidy$sra)

# save
readr::write_csv(runinfo_filt, str_replace(runinfo_path, "\\.runInfo\\.txt$", str_c(".", format(Sys.time(), "%Y%m%d"), ".filt", ".runInfo.txt")))
  