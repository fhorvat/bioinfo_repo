### INFO: 
### DATE: Sun Mar 04 11:21:12 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/pig_GV_Libechov_2018/Data/Raw/Links")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# mutate in rows for which cond is TRUE
mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
  
}

######################################################## READ DATA
# raw file path
raw_path <- "/common/RAW/Svoboda/pig_smallRNASeq_Libechov_2018"

# list sequenced files
seq_files <- list.files(path = raw_path, pattern = "*.fastq.gz", full.names = T)

######################################################## MAIN CODE
# clean names
names_df <- 
  tibble(path = seq_files) %>% 
  dplyr::mutate(ID = basename(path), 
                size = ifelse(str_detect(ID, "^L"), "large", "small"), 
                stage = "GV") %>% 
  dplyr::group_by(size) %>% 
  dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_out = str_c("s_", stage, "_", size, "_r", replicate, ".SE.txt.gz")) 
               
# get links script
links_df <- 
  names_df %>% 
  dplyr::select(path, sample_out) %>% 
  dplyr::mutate(sample_out = file.path(outpath, sample_out)) %>% 
  dplyr::mutate(link = str_c("ln -s", path, sample_out, sep = " ")) %$% 
  link %T>%
  readr::write_lines(., path = file.path(outpath, "make_links.sh"))
