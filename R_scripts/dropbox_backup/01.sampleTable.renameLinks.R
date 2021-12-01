### INFO:
### DATE: Sun Mar 04 11:21:12 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

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

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

# raw sequences path
fastq_path <- "/common/RAW/Svoboda/mESC_oocytes_2018"

######################################################## READ DATA
# read sample table
sample_df <- 
  readr::read_csv(file = file.path(outpath, "Libraries_sending_pools_180214_GeneCore_IDs_corrected.short.csv")) %>% 
  magrittr::set_colnames(., c("ID", "sample"))

######################################################## MAIN CODE
#### before mapping - create sample table
# get raw sequences table
fastq_df <-
  tibble(fastq_path = list.files(path = fastq_path, pattern = "*.txt.gz", recursive = T, full.names = T)) %>%
  mutate(fastq_ID = stringr::str_extract(basename(fastq_path),
                                         "MT[0-9]{1}|HET[0-9]{1}|SOM[0-9]{1}|WT[0-9]{1}|DXi[0-9]{1}|ICRF[0-9]{1}|ICRS[0-9]{1}|DXII.*sample[0-9]{1}") %>% 
           stringr::str_remove(., "_.*sample") %>% str_replace(., "DXII", "DXIIi"))

# clean sample table
sample_df_filt <-
  sample_df %>%
  dplyr::mutate(stage = rep(x = c("ESC", "GV", "MII", "ESC"), c(8, 16, 11, 8)),
                ID_clean = str_replace(ID, "ES ", "") %>% str_replace_all(., " ", "_") %>% str_replace(., "B6_WT", "WT_1")) %>%
  mutate_cond(stage == "ESC", sample_id = str_c(stage, ID_clean, str_replace_all(sample, " ", "_"), sep = "_")) %>%
  mutate_cond(stage == "GV", sample_id = str_replace_all(sample, " ", "_")) %>%
  mutate_cond(stage == "MII", sample_id = str_c(stage, "B6", ID_clean, sep = "_")) %>%
  dplyr::mutate(fastq_ID = str_replace_all(ID_clean, "_", ""),
                sample_id = str_c("s_", sample_id, ".SE"), 
                short_name = str_replace_all(sample_id, "DX_|i[0-9]{1}_|B6_|s_|.SE|MII_|ESC_|GV_|DXII_", ""), 
                genotype = str_replace(short_name, "_[0-9]{1}", ""), 
                sample_id = str_remove(sample_id, "(?<=DXII_i[0-9]{1}).*(?=.SE)")) %>%
  dplyr::right_join(., fastq_df, by = "fastq_ID") %>% 
  dplyr::select(ID, sample, sample_id, short_name, stage, genotype, fastq_path) %T>%
  readr::write_csv(., path = file.path(outpath, "mESC_oocytes_2018.sample_table.csv"))
  
# write script which creates renamed links
sample_df_filt %>%
  dplyr::mutate(sample_out = str_c(sample_id, ".txt.gz"), 
                link = str_c("ln -s ", fastq_path, " ./", sample_out)) %$%
  link %T>%
  readr::write_lines(., path = "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Raw/Links/rename_links.sh")


#### after mapping - add bam paths and library size to sample table
# mapped sequences path
mapped_path <- c("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Mapped/STAR_mm10", 
                 "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Mapped/STAR_mm10_20180601")

# read stats table
stats_df <- 
  purrr::map(list.files(path = mapped_path, pattern = "*stats_and_tracks.csv", full.names = T), readr::read_csv) %>% 
  dplyr::bind_rows(.)

# get mapped sequences table
bam_df <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".genome.Aligned.sortedByCoord.out.bam", ""))

# get library sizes
libSize_df <- 
  stats_df %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA)

# add to table
sample_df_filt <- 
  readr::read_csv(file = file.path(outpath, "mESC_oocytes_2018.sample_table.csv")) %>%
  dplyr::left_join(., bam_df, by = "sample_id") %>% 
  dplyr::left_join(., libSize_df, by = "sample_id") %>% 
  dplyr::select(ID:genotype, library_size, bam_path, fastq_path) %T>%
  readr::write_csv(., path = file.path(outpath, "mESC_oocytes_2018.sample_table.csv"))


