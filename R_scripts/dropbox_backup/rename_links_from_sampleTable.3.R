### INFO: 
### DATE: Sun Mar 04 11:21:12 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Hendrickson_2017_NatGenet_GSE72379/Data/Raw/Links")

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
# read sample table
sample_df <- 
  readr::read_csv(file = list.files(path = outpath, pattern = "*.csv")) %>% 
  dplyr::select(1:14) 

# list sequenced files
seq_files <- list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Raw/Links", pattern = "*.txt.gz")

######################################################## MAIN CODE
# clean sample table
sample_df_filt <- 
  sample_df %>% 
  dplyr::select(1:2) %>% 
  magrittr::set_colnames(., c("ID", "sample_name")) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::mutate(stage = rep(x = c("ESC", "GV", "MII"), c(8, 16, 10)), 
                ID = str_replace(ID, "ES ", ""), 
                ID = str_replace_all(ID, " ", "_")) %>% 
  mutate_cond(stage == "ESC", sample_out = str_c(stage, ID, sample_name, sep = "_")) %>%
  mutate_cond(stage == "GV", sample_out = str_replace_all(sample_name, " ", "_")) %>% 
  mutate_cond(stage == "MII", sample_out = str_c(stage, "B6", ID, sep = "_")) %>% 
  mutate_cond(ID == "B6_WT", sample_out = str_replace(sample_out, "_B6([^B6]*)$", "\\1")) %>% 
  dplyr::mutate(sample_out = str_c("s_", sample_out, ".SE.txt.gz"))
  
# sequence data.frame
seq_df <- 
  tibble(path = seq_files) %>% 
  dplyr::mutate(ID = str_replace_all(string = seq_files, ".*Malik_lane[0-9]|_sequence.txt.gz|ES", ""))

# join
rename_df <- 
  sample_df_filt %>% 
  dplyr::select(ID, sample_out) %>% 
  dplyr::mutate(ID = str_replace_all(ID, "_", "")) %>% 
  dplyr::right_join(seq_df, by = "ID") %>% 
  dplyr::mutate(mv = str_c("mv", path, sample_out, sep = " ")) %$% 
  mv %T>%
  readr::write_lines(., path = "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Raw/Links/02.rename_links.sh")
